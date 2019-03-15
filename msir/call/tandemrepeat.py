#!/usr/bin/env python

import logging
from multiprocessing import cpu_count
import os
import subprocess
from ..df.beddf import BedDataFrame
from ..df.samdf import SamDataFrame
from ..util.error import MsirError
from ..util.helper import fetch_abspath, fetch_executable, \
    validate_files_and_dirs


def extract_tandem_repeats_in_bams(bed_path, bam_paths, out_dir_path,
                                   index_bam=False, samtools=None,
                                   output_tsv=False, processes=None):
    validate_files_and_dirs(files=[bed_path, *bam_paths], dirs=[out_dir_path])
    bam_abspaths = [fetch_abspath(p) for p in bam_paths]
    samtools_path = samtools or fetch_executable('samtools')
    _validate_or_prepare_bam_indexes(
        bam_paths=bam_abspaths, index_bam=index_bam, samtools=samtools_path
    )
    beddf = BedDataFrame(path=fetch_abspath(bed_path))
    beddf.load()
    common_args = {
        'beddf': beddf, 'out_dir_path': fetch_abspath(out_dir_path),
        'csvsep': ('\t' if output_tsv else ','), 'samtools': samtools_path,
        'n_proc': int(processes or cpu_count())
    }
    for p in bam_abspaths:
        _extract_tandem_repeats(bam_path=p, **common_args)


def _validate_or_prepare_bam_indexes(bam_paths, index_bam=False,
                                     samtools=None):
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
            logger.debug('BAM/CRAM index files exist')
        elif index_bam:
            print('msir>\tIndex BAM/CRAM files.')
            for p in bam_paths_without_bai:
                args = [samtools, 'index', p]
                logger.debug('args: {}'.format(args))
                subprocess.run(args=args, check=True)
        else:
            raise MsirError(
                'BAM/CRAM indexes not found: {}'.format(
                    list(bam_paths_without_bai)
                )
            )


def _extract_tandem_repeats(bam_path, beddf, out_dir_path, samtools,
                            csvsep=',', n_proc=1):
    logger = logging.getLogger(__name__)
    print('msir>\tExtract repeats in:\t{}'.format(bam_path), flush=True)
    for id, bedline in beddf.df.iterrows():
        logger.debug('id: {}'.format(id))
        logger.debug('bedline:{0}{1}'.format(os.linesep, bedline))
        region = '{0}:{1}-{2}'.format(
            bedline['chrom'], bedline['chromStart'], bedline['chromEnd']
        )
        logger.debug('region: {}'.format(region))
        bamdf = SamDataFrame(path=bam_path, samtools=samtools, n_thread=n_proc)
        bamdf.load(add_view_args=[region])
        logger.debug('bamdf.df:{0}{1}'.format(os.linesep, bamdf.df))
