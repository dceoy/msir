#!/usr/bin/env python

from concurrent.futures import as_completed, ProcessPoolExecutor
from itertools import chain, product
import logging
import re
import numpy as np
import pandas as pd
from ..df.beddf import BedDataFrame
from ..util.helper import fetch_abspath, print_log, read_fasta, \
    validate_files_and_dirs


def identify_repeat_units_on_bed(bed_path, genome_fa_path, ru_tsv_path,
                                 max_unit_len=6, min_rep_times=3,
                                 ex_region_len=10, n_proc=8):
    logger = logging.getLogger(__name__)
    validate_files_and_dirs(files=[bed_path, genome_fa_path])
    print_log('Identify repeat units.')
    fa_seqs = read_fasta(fa_path=genome_fa_path)
    logger.debug('fa_seqs: {}'.format(fa_seqs))
    fa_seqs_len = {k: len(v.seq) for k, v in fa_seqs.items()}
    logger.debug('fa_seqs_len: {}'.format(fa_seqs_len))
    beddf = BedDataFrame(path=fetch_abspath(bed_path))
    beddf.load()
    beddf.df.assign(
        search_start=lambda d: d['chromStart'].apply(
            lambda i: max(0, i - ex_region_len)
        ),
        search_end=lambda d: d[['chrom', 'chromEnd']].apply(
            lambda r: min(r[0] - 1, r[1] + ex_region_len), axis=1
        )
    ).assign(
        search_seq=lambda d: d[['chrom', 'chromStart', 'chromEnd']].apply(
            lambda r: str(fa_seqs[r[0]].seq[r[1]:r[2]].upper())
        )
    )
    with ProcessPoolExecutor(max_workers=n_proc) as x:
        fs = [
            x.submit(
                _identify_repeat_unit, line, id, max_unit_len, min_rep_times
            ) for id, line in beddf.df.iterrows()
        ]
        df_u = pd.DataFrame([f.result() for f in as_completed(fs)], index='id')
    df_u.to_csv(fetch_abspath(ru_tsv_path), sep='\t')


def _identify_repeat_unit(region_dict, region_id, max_unit_len, min_rep_times):
    runits = chain.from_iterable([
        product('ATGC', repeat=(i + 1)) for i in range(max_unit_len)
    ])
    tr_list = [
        count_repeat_times(
            sequence=region_dict['search_seq'], repeat_unit=''.join(t),
            min_rep_times=min_rep_times, cut_end_len=0,
            start_pos=region_dict['chromStart']
        ) for t in runits[::-1]
    ]
    return {
        'id': region_id, **region_dict,
        **pd.DataFrame([d for d in tr_list if d]).pipe(
            lambda d: d[d['repeat_times'].idxmax()]
        ).to_dict()
    }


def count_repeat_times(sequence, repeat_unit, min_rep_times=1, cut_end_len=0,
                       start_pos=0):
    matches = [
        (*m.span(), m.group(0)) for m in re.finditer(
            r'({0}){{1},}'.format(repeat_unit, min_rep_times), sequence
        )
    ]
    if matches:
        return {
            'repeat_unit': repeat_unit, 'repeat_unit_size': len(repeat_unit),
            **pd.DataFrame(
                matches,
                columns=['repeat_start', 'repeat_end', 'repeat_sequence']
            ).pipe(
                lambda d: d[
                    (d['repeat_start'] > cut_end_len) &
                    (d['repeat_end'] < len(sequence) - cut_end_len)
                ]
            ).assign(
                repeat_times=lambda d: np.int32(
                    (d['repeat_end'] - d['repeat_start']) / len(repeat_unit)
                )
            ).pipe(
                lambda d: d.iloc[d['repeat_times'].idxmax()]
            ).assign(
                repeat_start_pos=lambda d: d['repeat_start'] + start_pos,
                repeat_end_pos=lambda d: d['repeat_end'] + start_pos
            ).drop(columns=['repeat_start', 'repeat_end'])
        }
    else:
        return dict()
