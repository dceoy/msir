#!/usr/bin/env python

from collections import OrderedDict
from concurrent.futures import as_completed, ProcessPoolExecutor
from itertools import chain, product
import logging
import os
from pprint import pformat
import re
import numpy as np
import pandas as pd
from ..df.beddf import BedDataFrame
from ..util.biotools import convert_bed_line_to_sam_region, read_fasta
from ..util.helper import fetch_abspath, print_log, validate_files_and_dirs


def identify_repeat_units_on_bed(bed_path, genome_fa_path, ru_tsv_path,
                                 max_unit_len=6, min_rep_times=3,
                                 ex_region_len=10, n_proc=8):
    logger = logging.getLogger(__name__)
    validate_files_and_dirs(files=[bed_path, genome_fa_path])
    print_log('Load FASTA data:\t{}'.format(genome_fa_path))
    fa_seqs = read_fasta(fa_path=genome_fa_path)
    fa_seqs_lens = {k: len(v.seq) for k, v in fa_seqs.items()}
    logger.debug('fa_seqs_lens:' + os.linesep + pformat(fa_seqs_lens))
    print_log('Load BED data:\t{}'.format(bed_path))
    beddf = BedDataFrame(path=fetch_abspath(bed_path))
    beddf.load()
    df_search = beddf.df.assign(
        search_start=lambda d: d['chromStart'].apply(
            lambda i: max(0, i - ex_region_len)
        ),
        search_end=lambda d: d[['chrom', 'chromEnd']].apply(
            lambda r: min(fa_seqs_lens[r[0]] - 1, r[1] + ex_region_len), axis=1
        )
    ).assign(
        search_seq=lambda d: d[['chrom', 'search_start', 'search_end']].apply(
            lambda r: str(fa_seqs[r[0]].seq[int(r[1]):int(r[2])].upper()),
            axis=1
        )
    )
    logger.debug('df_search:{0}{1}'.format(os.linesep, df_search))
    fa_seqs = None
    print_log('Identify repeat units on BED regions:')
    regex_dict = _make_str_regex_dict(
        max_unit_len=max_unit_len, min_rep_times=min_rep_times
    )
    logger.debug('len(regex_dict): {}'.format(len(regex_dict)))
    df_ru_list = []
    with ProcessPoolExecutor(max_workers=n_proc) as x:
        fs = [
            x.submit(_identify_repeat_unit, bedline, id, regex_dict)
            for id, bedline in df_search.iterrows()
        ]
        try:
            for f in as_completed(fs):
                res = f.result()
                _print_state_line(res=res)
                df_ru_list.append(res['df'])
        except Exception as e:
            [f.cancel() for f in fs]
            x.shutdown(wait=False)
            raise e
        else:
            df_ru = pd.concat(
                df_ru_list
            ).reset_index().set_index('id').sort_index()
            logger.debug('df_ru:{0}{1}'.format(os.linesep, df_ru))
    ru_tsv_abspath = fetch_abspath(ru_tsv_path)
    print_log('Write repeat units data:\t{}'.format(ru_tsv_abspath))
    df_ru.to_csv(ru_tsv_abspath, index=False, sep='\t')


def _print_state_line(res):
    d = res['df'].iloc[0].to_dict() if res['df'].size else dict()
    line = '  {0:<25}\t{1:<10}\t{2}'.format(
        res['region'],
        ('{0}x{1}'.format(d['repeat_times'], d['repeat_unit']) if d else '-'),
        (d['repeat_seq'] if d else '-')
    )
    print(line, flush=True)


def _make_str_regex_dict(max_unit_len=6, min_rep_times=1, bases='ACGT'):
    units = [
        ''.join(t) for t in chain.from_iterable([
            product(bases, repeat=i) for i in range(max_unit_len, 0, -1)
        ])
    ]
    return OrderedDict([
        (s, compile_str_regex(s, min_rep_times)) for s in units
    ])


def compile_str_regex(repeat_unit, min_rep_times=1):
    return re.compile(r'(%s){%d,}' % (repeat_unit, min_rep_times))


def _identify_repeat_unit(region_dict, region_id, regex_dict):
    df_lr = extract_longest_repeat_df(
        sequence=region_dict['search_seq'], regex_dict=regex_dict,
        cut_end_len=0, start_pos=region_dict['chromStart']
    )
    return {
        'region': convert_bed_line_to_sam_region(region_dict),
        'df': (
            df_lr.assign(id=region_id, **region_dict,).set_index([
                'id', *region_dict.keys()
            ]) if df_lr.size else pd.DataFrame()
        )
    }


def extract_longest_repeat_df(sequence, regex_dict, cut_end_len=0,
                              start_pos=0):
    hits = chain.from_iterable([
        [(*m.span(), m.group(0), u) for m in r.finditer(sequence)]
        for u, r in regex_dict.items() if u in sequence
    ])
    if hits:
        return pd.DataFrame(
            hits, columns=['start_idx', 'end_idx', 'repeat_seq', 'repeat_unit']
        ).assign(
            repeat_seq_size=lambda d: d['end_idx'] - d['start_idx'],
            repeat_unit_size=lambda d: d['repeat_unit'].apply(len),
            repeat_start=lambda d: d['start_idx'] + start_pos,
            repeat_end=lambda d: d['end_idx'] + start_pos
        ).assign(
            repeat_times=lambda d:
            np.int32(d['repeat_seq_size'] / d['repeat_unit_size'])
        ).pipe(
            lambda d: d.iloc[[d['repeat_seq_size'].idxmax()]]
        ).pipe(
            lambda d: d[
                (d['start_idx'] >= cut_end_len) &
                (d['end_idx'] <= len(sequence) - cut_end_len)
            ].drop(columns=['start_idx', 'end_idx'])
        )
    else:
        return pd.DataFrame()
