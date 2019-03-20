#!/usr/bin/env python

from collections import OrderedDict
from concurrent.futures import as_completed, ProcessPoolExecutor
from itertools import chain
import logging
import os
from pprint import pformat
import re
import pandas as pd
from ..util.biotools import convert_bed_line_to_sam_region, \
    iterate_unique_repeat_units, read_fasta, read_bed
from ..util.helper import fetch_abspath, print_log, validate_files_and_dirs


def identify_repeat_units_on_bed(bed_path, genome_fa_path, ru_tsv_path,
                                 max_unit_len=6, min_rep_times=3,
                                 ex_region_len=10, n_proc=8):
    validate_files_and_dirs(files=[bed_path, genome_fa_path])
    print_log('Load input data:')
    df_exbed = _make_extented_bed_df(
        bed_path=bed_path, genome_fa_path=genome_fa_path,
        ex_region_len=ex_region_len
    )
    print_log('Compile regular expression patterns:', end='')
    regex_patterns = _compile_repeat_unit_regex_patterns(
        max_unit_len=max_unit_len, min_rep_times=min_rep_times
    )
    print('\t{}'.format(len(regex_patterns)), flush=True)
    print_log('Identify repeat units on BED regions:')
    df_ru = _make_repeat_unit_df(
        df_exbed=df_exbed, regex_patterns=regex_patterns, n_proc=n_proc
    )
    print_log('Write repeat units data:\t{}'.format(ru_tsv_path))
    df_ru.to_csv(fetch_abspath(ru_tsv_path), index=False, sep='\t')


def _make_extented_bed_df(bed_path, genome_fa_path, ex_region_len=10):
    logger = logging.getLogger(__name__)
    print('  FASTA:\t{}'.format(genome_fa_path), flush=True)
    ref_genome = read_fasta(path=genome_fa_path)
    seq_lengths = {k: len(v.seq) for k, v in ref_genome.items()}
    logger.debug('seq_lengths:' + os.linesep + pformat(seq_lengths))
    print('  BED:\t{}'.format(bed_path), flush=True)
    df_exbed = read_bed(path=bed_path).assign(
        search_start=lambda d: d['chromStart'].apply(
            lambda i: max(0, i - ex_region_len)
        ),
        search_end=lambda d: d[['chrom', 'chromEnd']].apply(
            lambda r: min(seq_lengths[r[0]] - 1, r[1] + ex_region_len), axis=1
        )
    ).assign(
        search_seq=lambda d: d[['chrom', 'search_start', 'search_end']].apply(
            lambda r: str(ref_genome[r[0]].seq[int(r[1]):int(r[2])].upper()),
            axis=1
        )
    )
    logger.debug('df_exbed:{0}{1}'.format(os.linesep, df_exbed))
    return df_exbed


def _compile_repeat_unit_regex_patterns(max_unit_len=6, min_rep_times=1):
    return OrderedDict([
        (s, compile_str_regex(s, min_rep_times=min_rep_times))
        for s in iterate_unique_repeat_units(max_unit_len=max_unit_len)
    ])


def compile_str_regex(repeat_unit, min_rep_times=1):
    return re.compile(r'(%s){%d,}' % (repeat_unit, min_rep_times))


def _make_repeat_unit_df(df_exbed, regex_patterns, n_proc=8):
    ppx = ProcessPoolExecutor(max_workers=n_proc)
    fs = [
        ppx.submit(_identify_repeat_unit, bedline, id, regex_patterns)
        for id, bedline in df_exbed.iterrows()
    ]
    df_ru = pd.DataFrame()
    try:
        for f in as_completed(fs):
            res = f.result()
            _print_state_line(res=res)
            df_ru = df_ru.append(res['df'])
    except Exception as e:
        [f.cancel() for f in fs]
        ppx.shutdown(wait=False)
        raise e
    else:
        ppx.shutdown(wait=True)
        df_ru = df_ru.set_index('id').sort_index()
    finally:
        print(df_ru, flush=True)
    return df_ru


def _identify_repeat_unit(bedline, region_id, regex_patterns):
    return {
        'region': convert_bed_line_to_sam_region(bedline),
        'df': extract_longest_repeat_df(
            sequence=bedline['search_seq'], regex_patterns=regex_patterns,
            cut_end_len=0, start_pos=bedline['chromStart']
        ).assign(
            id=region_id, **bedline
        ).set_index(['id', *bedline.keys()]).reset_index()
    }


def extract_longest_repeat_df(sequence, regex_patterns, cut_end_len=0,
                              start_pos=0):
    hits = chain.from_iterable([
        [(*m.span(), m.group(0), u) for m in r.finditer(sequence)]
        for u, r in regex_patterns.items() if u in sequence
    ])
    if hits:
        return pd.DataFrame(
            hits, columns=['start_idx', 'end_idx', 'repeat_seq', 'repeat_unit']
        ).assign(
            repeat_seq_length=lambda d: d['end_idx'] - d['start_idx'],
            repeat_unit_length=lambda d: d['repeat_unit'].apply(len),
            repeat_start=lambda d: d['start_idx'] + start_pos,
            repeat_end=lambda d: d['end_idx'] + start_pos
        ).assign(
            repeat_times=lambda d:
            (d['repeat_seq_length'] / d['repeat_unit_length']).astype(int)
        ).pipe(
            lambda d: d[d['repeat_seq_length'] == d['repeat_seq_length'].max()]
        ).pipe(
            lambda d: d[
                (d['start_idx'] >= cut_end_len) &
                (d['end_idx'] <= len(sequence) - cut_end_len)
            ].drop(columns=['start_idx', 'end_idx'])
        ).reset_index().pipe(
            lambda d: d.iloc[[d['repeat_times'].idxmax()]]
        )
    else:
        return pd.DataFrame()


def _print_state_line(res):
    d = res['df'].iloc[0].to_dict() if res['df'].size else dict()
    line = '  {0:<25}\t{1:<10}\t{2}'.format(
        res['region'],
        ('{0}x{1}'.format(d['repeat_times'], d['repeat_unit']) if d else '-'),
        (d['repeat_seq'] if d else '-')
    )
    print(line, flush=True)
