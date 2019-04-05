#!/usr/bin/env python

from collections import OrderedDict
from concurrent.futures import as_completed, ProcessPoolExecutor
from itertools import chain
import logging
import os
import traceback
from pprint import pformat
import re
import pandas as pd
from ..util.biotools import convert_bed_line_to_sam_region, \
    iterate_unique_repeat_units, read_fasta, read_bed
from ..util.helper import fetch_abspath, print_log, validate_files_and_dirs


def identify_repeat_units_on_bed(bed_path, genome_fa_path, trunit_tsv_path,
                                 max_unit_len=6, min_rep_times=3,
                                 min_rep_len=10, flanking_len=10,
                                 ex_region_len=20, n_proc=8):
    validate_files_and_dirs(files=[bed_path, genome_fa_path])
    print_log('Load input data:')
    df_exbed = _make_extented_bed_df(
        bed_path=bed_path, genome_fa_path=genome_fa_path,
        ex_region_len=int(ex_region_len)
    )
    print_log('Compile regular expression patterns:', end='')
    regex_patterns = _compile_repeat_unit_regex_patterns(
        max_unit_len=int(max_unit_len), min_rep_times=int(min_rep_times)
    )
    print('\t{}'.format(len(regex_patterns['patterns'])), flush=True)
    print_log('Identify repeat units on BED regions:')
    df_ru = _make_repeat_unit_df(
        df_exbed=df_exbed, regex_patterns=regex_patterns,
        min_rep_len=int(min_rep_len), flanking_len=int(flanking_len),
        n_proc=n_proc
    )
    if df_ru.size:
        print_log('Write repeat units data:\t{}'.format(trunit_tsv_path))
        df_ru.to_csv(fetch_abspath(trunit_tsv_path), sep='\t')
    else:
        print_log('Failed to identify repeat units.')


def _make_extented_bed_df(bed_path, genome_fa_path, ex_region_len=10):
    logger = logging.getLogger(__name__)
    print('  FASTA:\t{}'.format(genome_fa_path), flush=True)
    ref_genome = read_fasta(path=genome_fa_path)
    seq_lens = {k: len(v.seq) for k, v in ref_genome.items()}
    logger.debug('seq_lens:' + os.linesep + pformat(seq_lens))
    print('  BED:\t{}'.format(bed_path), flush=True)
    df_exbed = read_bed(
        path=bed_path
    )[['chrom', 'chromStart', 'chromEnd']].assign(
        search_start=lambda d: d['chromStart'].apply(
            lambda i: max(0, i - ex_region_len)
        ),
        search_end=lambda d: d[['chrom', 'chromEnd']].apply(
            lambda r: min(seq_lens[r[0]] - 1, r[1] + ex_region_len), axis=1
        )
    ).assign(
        search_seq=lambda d: d[['chrom', 'search_start', 'search_end']].apply(
            lambda r: str(ref_genome[r[0]].seq[int(r[1]):int(r[2])].upper()),
            axis=1
        )
    )
    logger.debug('df_exbed:{0}{1}'.format(os.linesep, df_exbed))
    return df_exbed


def _compile_repeat_unit_regex_patterns(max_unit_len=6, **kwargs):
    return {
        'patterns': OrderedDict([
            (s, compile_str_regex(repeat_unit=s, **kwargs))
            for s in iterate_unique_repeat_units(max_unit_len=max_unit_len)
        ]),
        'max_unit_len': max_unit_len,
        **kwargs
    }


def compile_str_regex(repeat_unit, min_rep_times=1, left_seq='', right_seq=''):
    if left_seq or right_seq:
        return re.compile(
            r'%s(%s){%d,}%s' %
            (left_seq, repeat_unit, min_rep_times, right_seq)
        )
    else:
        return re.compile(r'(%s){%d,}' % (repeat_unit, min_rep_times))


def _make_repeat_unit_df(df_exbed, regex_patterns, min_rep_len=10,
                         flanking_len=0, n_proc=8):
    logger = logging.getLogger(__name__)
    ppx = ProcessPoolExecutor(max_workers=n_proc)
    fs = [
        ppx.submit(
            _identify_repeat_unit, bedline, id, regex_patterns, min_rep_len,
            flanking_len
        ) for id, bedline in df_exbed.iterrows()
    ]
    try:
        df_ru = pd.concat(
            [f.result() for f in as_completed(fs)], sort=False
        )
    except Exception as e:
        logger.error(os.linesep + traceback.format_exc())
        ppx.shutdown(wait=False)
        raise e
    else:
        logger.debug('df_ru:{0}{1}'.format(os.linesep, df_ru))
        ppx.shutdown(wait=True)
    if df_ru.size:
        return df_ru.sort_values(by='bed_id').drop(columns='bed_id')
    else:
        return df_ru


def _identify_repeat_unit(bedline, bed_id, regex_patterns, min_rep_len,
                          flanking_len):
    seq = bedline['search_seq']
    df_line = extract_longest_repeat_df(
        sequence=seq, regex_patterns=regex_patterns, min_rep_len=min_rep_len,
        flanking_len=flanking_len, start_pos=bedline['chromStart']
    ).pipe(
        lambda d: (
            d.assign(
                bed_id=bed_id, **bedline,
                left_seq=lambda d: d[['start_x', 'end_x']].apply(
                    lambda r: seq[(r[0] - flanking_len):r[0]], axis=1
                ),
                right_seq=lambda d: d[['start_x', 'end_x']].apply(
                    lambda r: seq[r[1]:(r[1] + flanking_len)], axis=1
                )
            ).set_index(['chrom', 'chromStart', 'chromEnd'])[[
                'bed_id', 'repeat_start', 'repeat_end', 'repeat_unit',
                'repeat_unit_length', 'repeat_times', 'repeat_seq_length',
                'left_seq', 'repeat_seq', 'right_seq', 'search_start',
                'search_end', 'search_seq'
            ]].sort_index()
            if d.size else d
        )
    )
    _print_state_line(
        region=convert_bed_line_to_sam_region(bedline), df=df_line
    )
    return df_line


def extract_longest_repeat_df(sequence, regex_patterns, min_rep_len=2,
                              flanking_len=0, start_pos=0):
    lsl = len(regex_patterns.get('left_seq') or '')
    rsl = len(regex_patterns.get('right_seq') or '')
    hits = [
        t for t in chain.from_iterable([
            [(m.group(0), u, *m.span()) for m in r.finditer(sequence)]
            for u, r in regex_patterns['patterns'].items() if u in sequence
        ])
        if (len(t[0]) >= sum([min_rep_len, lsl, rsl]) and
            t[2] >= flanking_len and t[3] + flanking_len <= len(sequence))
    ]
    if not hits:
        return pd.DataFrame()
    else:
        return pd.DataFrame(
            hits, columns=['repeat_seq', 'repeat_unit', 'start_x', 'end_x']
        ).assign(
            repeat_start=lambda d: d['start_x'] + lsl + start_pos,
            repeat_end=lambda d: d['end_x'] - rsl + start_pos,
            repeat_unit_length=lambda d: d['repeat_unit'].apply(len),
            left_seq=lambda d: (d['repeat_seq'].str[:lsl] if lsl else ''),
            right_seq=lambda d: (d['repeat_seq'].str[-rsl:] if rsl else '')
        ).assign(
            repeat_seq_length=lambda d: d['repeat_end'] - d['repeat_start']
        ).assign(
            repeat_times=lambda d:
            (d['repeat_seq_length'] / d['repeat_unit_length']).astype(int)
        ).pipe(
            lambda d: d[d['repeat_seq_length'] == d['repeat_seq_length'].max()]
        ).reset_index(drop=True).pipe(
            lambda d: d.iloc[[d['repeat_times'].idxmax()]]
        )


def _print_state_line(region, df):
    d = df.iloc[0].to_dict() if df.size else dict()
    line = '  {0:<25}\t{1:<10}\t{2}'.format(
        region,
        ('{0}x{1}'.format(d['repeat_times'], d['repeat_unit']) if d else '-'),
        (d['left_seq'] + d['repeat_seq'] + d['right_seq'] if d else '-')
    )
    print(line, flush=True)
