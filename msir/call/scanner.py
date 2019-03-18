#!/usr/bin/env python

from concurrent.futures import as_completed, ProcessPoolExecutor
import logging
import os
import pandas as pd
from ..util.helper import fetch_abspath, print_log, validate_files_and_dirs
from ..util.biotools import fetch_bed_region_str, \
    validate_or_prepare_bam_indexes, view_bam_lines
from .identifier import compile_str_regex, extract_longest_repeat_df


def scan_tandem_repeats_in_reads(bam_paths, ru_tsv_path, out_dir_path,
                                 out_file_path=None, index_bam=False,
                                 samtools=None, cut_end_len=10,
                                 use_csv_format=False, n_proc=8):
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
    regex_dict = {
        s: compile_str_regex(repeat_unit=s, min_rep_times=1)
        for s in set(df_ru['repeat_unit'])
    }
    logger.debug('list(regex_dict.keys()): {}'.format(list(regex_dict.keys())))
    tabfmt = (
        {'ext': 'csv', 'sep': ','} if use_csv_format
        else {'ext': 'tsv', 'sep': '\t'}
    )
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
                _print_state_line(res=res, bam_path=p)
                df_read_list.append(res['df'])
        df_read = pd.concat(df_read_list).set_index('id').sort_index()
        logger.debug('df_read:{0}{1}'.format(os.linesep, df_read))
        res_tsv_path = fetch_abspath(
            out_file_path or os.path.join(
                out_dir_path,
                '{0}.repears.{1}'.format(os.path.basename(p), tabfmt['ext'])
            )
        )
        print_log('Write results into:\t{}'.format(res_tsv_path))
        if out_file_path:
            df_read.pipe(
                lambda d: d.assign(data_path=p)[['data_path', *d.columns]]
            ).to_csv(
                res_tsv_path, index=False, sep=tabfmt['sep'], mode='a',
                header=(not os.path.isfile(res_tsv_path))
            )
        else:
            df_read.to_csv(res_tsv_path, index=False, sep=tabfmt['sep'])
    logger.debug('all the processes done.')


def _print_state_line(res, bam_path):
    if res['df'].size:
        for i, r in res['df'].iterrows():
            line = '  {0}\t{1:<25}\t{2:<10}\t{3}'.format(
                os.path.basename(bam_path), res['region'],
                '{0}x{1}'.format(r['read_repeat_times'], r['repeat_unit']),
                r['read_repeat_times_count']
            )
            print(line, flush=True)


def _count_repeats_in_reads(bam_path, tsvline, id, regex_dict, cut_end_len,
                            samtools):
    region = fetch_bed_region_str(**tsvline)
    bed_cols = ['chrom', 'chromStart', 'chromEnd']
    ru = tsvline['repeat_unit']
    df_bam = pd.concat([
        extract_longest_repeat_df(
            sequence=d['SEQ'], regex_dict={ru: regex_dict[ru]},
            cut_end_len=cut_end_len, start_pos=d['POS']
        ) for d in view_bam_lines(
            bam_path=bam_path, regions=[region], samtools_path=samtools
        )
    ])
    df_rt = (
        df_bam['repeat_times'].value_counts().to_frame(
            name='read_repeat_times_count'
        ).reset_index().rename(
            columns={'index': 'read_repeat_times'}
        ).assign(
            id=id, **{k: tsvline[k] for k in bed_cols}, repeat_unit=ru
        ).sort_values(
            'read_repeat_times_count', ascending=False
        )[[
            'id', *bed_cols, 'repeat_unit', 'read_repeat_times',
            'read_repeat_times_count'
        ]]
        if df_bam.size else pd.DataFrame()
    )
    return {'region': region, 'df': df_rt}
