#!/usr/bin/env python

from concurrent.futures import as_completed, ProcessPoolExecutor
import logging
import os
import subprocess
import pandas as pd
from ..util.error import MsirError
from ..util.helper import fetch_abspath, fetch_executable, print_log, \
    run_and_parse_subprocess, validate_files_and_dirs


def scan_tandem_repeats_in_reads(bam_paths, ru_tsv_path, out_dir_path,
                                 index_bam=False, samtools=None,
                                 cut_end_len=10, output_csv=False, n_proc=8):
    validate_files_and_dirs(
        files=[ru_tsv_path, *bam_paths], dirs=[out_dir_path]
    )
    bam_abspaths = [fetch_abspath(p) for p in bam_paths]
    samtools_path = samtools or fetch_executable('samtools')
    _validate_or_prepare_bam_indexes(
        bam_paths=bam_abspaths, index_bam=index_bam, n_proc=n_proc,
        samtools=samtools_path
    )
    df_ru = pd.read_csv(ru_tsv_path, sep='\t')
    out_dir_abspath = fetch_abspath(out_dir_path)
    table_ext = 'csv' if output_csv else 'tsv'
    for p in bam_abspaths:
        print_log('Extract repeats in: {}'.format(p))
        with ProcessPoolExecutor(max_workers=n_proc) as x:
            fs = [
                x.submit(_count_repeats_in_reads, p, line, id, samtools_path)
                for id, line in df_ru.iterrows()
            ]
            df_res = pd.DataFrame(
                [f.result() for f in as_completed(fs)], index='id'
            )
        result_tsv_path = os.path.join(
            out_dir_abspath, '{0}.msir.{1}'.format(p, table_ext)
        )
        df_res.to_csv(
            result_tsv_path, mode='a',
            header=(not os.path.isfile(result_tsv_path)),
            sep={'csv': ',', 'tsv': '\t'}[table_ext]
        )


def _validate_or_prepare_bam_indexes(bam_paths, index_bam, n_proc, samtools):
    logger = logging.getLogger(__name__)
    invalid_exts = [p for p in bam_paths if not p.endswith(('.bam', '.cram'))]
    if invalid_exts:
        raise MsirError('invalid BAM/CRAM paths: {}'.format(invalid_exts))
    else:
        bai_paths = {
            p: (p + '.bai' if p.endswith('.bam') else p + '.crai')
            for p in bam_paths
        }
        bam_paths_without_bai = [
            k for k, v in bai_paths.items() if not os.path.isfile(v)
        ]
        if not bam_paths_without_bai:
            pass
        elif index_bam:
            print_log('Index BAM/CRAM files.')
            for p in bam_paths_without_bai:
                args = [samtools, 'index', '-@', n_proc, p]
                logger.debug('args: {}'.format(args))
                subprocess.run(args=args, check=True)
        else:
            raise MsirError(
                'BAM/CRAM indexes not found: {}'.format(
                    list(bam_paths_without_bai)
                )
            )


def _count_repeats_in_reads(bam_path, tsvline, id, samtools):
    logger = logging.getLogger(__name__)
    logger.debug('id: {}'.format(id))
    logger.debug('tsvline:{0}{1}'.format(os.linesep, tsvline))
    region = '{0}:{1}-{2}'.format(
        tsvline['chrom'], tsvline['chromStart'], tsvline['chromEnd']
    )
    logger.debug('region: {}'.format(region))
    args = [samtools, 'view', '-F', '4', bam_path, region]
    logger.debug('args: {}'.format(args))
    for r in run_and_parse_subprocess(args=args):
        logger.debug(r)
    return dict()
