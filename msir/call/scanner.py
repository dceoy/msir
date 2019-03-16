#!/usr/bin/env python

from concurrent.futures import as_completed, ProcessPoolExecutor
import logging
import os
from pprint import pformat
import pandas as pd
from ..util.helper import fetch_abspath, fetch_bed_region_str, print_log, \
    validate_files_and_dirs
from ..util.samtools import validate_or_prepare_bam_indexes, view_bam_lines
from .identifier import count_repeat_times


def scan_tandem_repeats_in_reads(bam_paths, ru_tsv_path, out_dir_path,
                                 index_bam=False, samtools=None,
                                 cut_end_len=10, output_csv=False, n_proc=8):
    validate_files_and_dirs(
        files=[ru_tsv_path, *bam_paths], dirs=[out_dir_path]
    )
    bam_abspaths = [fetch_abspath(p) for p in bam_paths]
    validate_or_prepare_bam_indexes(
        bam_paths=bam_abspaths, index_bam=index_bam, n_proc=n_proc,
        samtools_path=samtools
    )
    print_log('Load repeat units data:\t{}'.format(ru_tsv_path))
    df_ru = pd.read_csv(ru_tsv_path, sep='\t')
    out_dir_abspath = fetch_abspath(out_dir_path)
    table_ext = 'csv' if output_csv else 'tsv'
    for p in bam_abspaths:
        print_log('Extract tandem repeats on reads:\t{}'.format(p))
        with ProcessPoolExecutor(max_workers=n_proc) as x:
            fs = [
                x.submit(
                    _count_repeats_in_reads, p, line, id, samtools
                ) for id, line in df_ru.iterrows()
            ]
            df_res = pd.concat(
                [f.result() for f in as_completed(fs)], index='id'
            )
        result_tsv_path = os.path.join(
            out_dir_abspath, '{0}.msir.{1}'.format(p, table_ext)
        )
        print_log('Write results into:\t{}'.format(result_tsv_path))
        df_res.to_csv(
            result_tsv_path, mode='a',
            header=(not os.path.isfile(result_tsv_path)),
            sep={'csv': ',', 'tsv': '\t'}[table_ext]
        )


def _count_repeats_in_reads(bam_path, tsvline, id, cut_end_len,
                            samtools):
    logger = logging.getLogger(__name__)
    logger.debug('tsvline:{0}{1}'.format(os.linesep, tsvline))
    repeats = [
        count_repeat_times(
            sequence=d['SEQ'], repeat_unit=tsvline['repeat_unit'],
            min_rep_times=1, cut_end_len=cut_end_len, start_pos=d['POS']
        ) for d in view_bam_lines(
            bam_path=bam_path, regions=[fetch_bed_region_str(**tsvline)],
            options=['-F', '4'], samtools_path=samtools
        )
    ]
    logger.debug('repeats:' + os.linesep + pformat(repeats))
    return pd.Series([
        d['repeat_times'] for d in repeats
    ]).value_counts().to_frame(name='count').reset_index().rename(
        columns={'index': 'repeat_times_in_read'}
    ).assign(id=id, **tsvline)
