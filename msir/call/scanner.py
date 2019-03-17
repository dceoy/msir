#!/usr/bin/env python

from collections import OrderedDict
from concurrent.futures import as_completed, ProcessPoolExecutor
import logging
import os
from pprint import pformat
import pandas as pd
from ..util.helper import fetch_abspath, fetch_bed_region_str, print_log, \
    validate_files_and_dirs
from ..util.samtools import validate_or_prepare_bam_indexes, view_bam_lines
from .identifier import compile_str_regex, extract_longest_repeat_df


def scan_tandem_repeats_in_reads(bam_paths, ru_tsv_path, out_dir_path,
                                 index_bam=False, samtools=None,
                                 cut_end_len=10, output_csv=False, n_proc=8):
    logger = logging.getLogger(__name__)
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
    regex_dict = OrderedDict([
        compile_str_regex(repeat_unit=s, min_rep_times=1)
        for s in set(df_ru['repeat_unit'])
    ])
    logger.debug('regex_dict:' + os.linesep + pformat(regex_dict))
    out_dir_abspath = fetch_abspath(out_dir_path)
    table_ext = 'csv' if output_csv else 'tsv'
    for p in bam_abspaths:
        print_log('Extract tandem repeats on reads:\t{}'.format(p))
        df_read_list = []
        with ProcessPoolExecutor(max_workers=n_proc) as x:
            fs = [
                x.submit(
                    _count_repeats_in_reads, p, line, id, regex_dict,
                    cut_end_len, samtools
                ) for id, line in df_ru.iterrows()
            ]
            for f in as_completed(fs):
                res = f.result()
                d = res['df'].assign(
                    region=res['region']
                )[['region', 'repeat_times_in_read', 'repeat_times_count']]
                print(d.to_string(index=False, header=False), flush=True)
                df_read_list.append(res['df'])
        df_read = pd.concat(df_read_list).set_index('id').sort_index()
        logger.debug('df_read:{0}{1}'.format(os.linesep, df_read))
        res_tsv_path = os.path.join(
            out_dir_abspath, '{0}.msir.{1}'.format(p, table_ext)
        )
        print_log('Write results into:\t{}'.format(res_tsv_path))
        df_read.to_csv(
            res_tsv_path, index=False, sep={'csv': ',', 'tsv': '\t'}[table_ext]
        )


def _count_repeats_in_reads(bam_path, tsvline, id, regex_dict, cut_end_len,
                            samtools):
    region = fetch_bed_region_str(**tsvline)
    bed_cols = ['chrom', 'chromStart', 'chromEnd']
    return {
        'region': region,
        'df': pd.concat([
            extract_longest_repeat_df(
                sequence=d['SEQ'],
                regex_dict={tsvline['repeat_unit']: regex_dict['repeat_unit']},
                cut_end_len=cut_end_len, start_pos=d['POS']
            )['repeat_times']
            for d in view_bam_lines(
                bam_path=bam_path, regions=[region], options=['-F', '4'],
                samtools_path=samtools
            )
        ]).value_counts().to_frame(
            name='repeat_times_count'
        ).reset_index().rename(
            columns={'index': 'repeat_times_in_read'}
        ).assign(
            id=id, **{k: tsvline[k] for k in bed_cols}
        ).sort_values(
            'repeat_times_count', ascending=False
        )[['id', *bed_cols, 'repeat_times_in_read', 'repeat_times_count']]
    }
