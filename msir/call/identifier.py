#!/usr/bin/env python

from concurrent.futures import as_completed, ProcessPoolExecutor
from itertools import chain, product
import logging
import os
from pprint import pformat
import re
import numpy as np
import pandas as pd
from ..df.beddf import BedDataFrame
from ..util.helper import fetch_abspath, fetch_bed_region_str, print_log, \
    read_fasta, validate_files_and_dirs


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
        search_seq=lambda d: d[['chrom', 'chromStart', 'chromEnd']].apply(
            lambda r: str(fa_seqs[r[0]].seq[int(r[1]):int(r[2])].upper()),
            axis=1
        )
    )
    logger.debug('df_search:{0}{1}'.format(os.linesep, df_search))
    fa_seqs = None
    df_ru = pd.DataFrame()
    print_log('Identify repeat units on BED regions:')
    with ProcessPoolExecutor(max_workers=n_proc) as x:
        fs = [
            x.submit(
                _identify_repeat_unit, line, id, max_unit_len, min_rep_times
            ) for id, line in df_search.iterrows()
        ]
        for f in as_completed(fs):
            res = f.result()
            rs = (res['df'].iloc[0]['repeat_seq'] if res['df'].size else '-')
            print('{0}\t>>\t{1}'.format(res['region'], rs), flush=True)
            df_ru.append(res['df'])
    df_ru.set_index('id').pipe(
        lambda d: d[[
            *df_search.columns,
            *[c for c in d.columns if c not in df_search.columns]
        ]]
    ).sort_index()
    logger.debug('df_ru:{0}{1}'.format(os.linesep, df_ru))
    ru_tsv_abspath = fetch_abspath(ru_tsv_path)
    print_log('Write repeat units data:\t{}'.format(ru_tsv_abspath))
    df_ru.to_csv(ru_tsv_abspath, index=False, sep='\t')


def _identify_repeat_unit(region_dict, region_id, max_unit_len, min_rep_times):
    ru_list = [
        count_repeat_times(
            sequence=region_dict['search_seq'], repeat_unit=''.join(t),
            min_rep_times=min_rep_times, cut_end_len=0,
            start_pos=region_dict['chromStart']
        ) for t in chain.from_iterable([
            product('ATGC', repeat=(i + 1)) for i in range(max_unit_len)
        ])
    ]
    ru_list_filtered = [d for d in ru_list[::-1] if d]
    return {
        'region': fetch_bed_region_str(**region_dict),
        'df': (
            pd.DataFrame(ru_list_filtered).pipe(
                lambda d: d.iloc[[d['repeat_times'].idxmax()]]
            ).assign(
                id=region_id, **region_dict,
            ) if ru_list_filtered else pd.DataFrame()
        )
    }


def count_repeat_times(sequence, repeat_unit, min_rep_times=1, cut_end_len=0,
                       start_pos=0):
    matches = [
        (*m.span(), m.group(0)) for m in re.finditer(
            (r'(%s){%d,}' % (repeat_unit, min_rep_times)), sequence
        )
    ]
    if matches:
        df_m = pd.DataFrame(
            matches, columns=['start_idx', 'end_idx', 'repeat_seq']
        ).pipe(
            lambda d: d[
                (d['start_idx'] >= cut_end_len) &
                (d['end_idx'] <= len(sequence) - cut_end_len)
            ]
        )
        if df_m.shape[0]:
            return {
                'repeat_unit': repeat_unit,
                'repeat_unit_size': len(repeat_unit),
                **df_m.assign(
                    repeat_times=lambda d: np.int32(
                        (d['end_idx'] - d['start_idx']) / len(repeat_unit)
                    ),
                    repeat_start=lambda d: d['start_idx'] + start_pos,
                    repeat_end=lambda d: d['end_idx'] + start_pos - 1
                ).drop(
                    columns=['start_idx', 'end_idx']
                ).pipe(
                    lambda d: d.iloc[d['repeat_times'].idxmax()]
                ).to_dict()
            }
        else:
            return dict()
    else:
        return dict()
