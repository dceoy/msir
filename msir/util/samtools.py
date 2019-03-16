#!/usr/bin/env python

import io
import logging
import os
import subprocess
import pandas as pd
from ..util.error import MsirError
from ..util.helper import fetch_executable, print_log, run_and_parse_subprocess


def validate_or_prepare_bam_indexes(bam_paths, index_bam=False, n_proc=8,
                                    samtools_path=None):
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
            samtools = samtools_path or fetch_executable('samtools')
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


def view_bam_lines(bam_path, regions=[], options=[], samtools_path=None):
    logger = logging.getLogger(__name__)
    samtools = samtools_path or fetch_executable('samtools')
    args = [samtools, 'view', *options, bam_path, *regions]
    logger.debug('args: {}'.format(args))
    for r in run_and_parse_subprocess(args=args):
        yield _parse_sam_line(line=r)


def _parse_sam_line(line):
    fixed_cols = [
        'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT',
        'TLEN', 'SEQ', 'QUAL'
    ]
    fixed_col_dtypes = {
        'QNAME': str, 'FLAG': int, 'RNAME': str, 'POS': int, 'MAPQ': int,
        'CIGAR': str, 'RNEXT': str, 'PNEXT': int, 'TLEN': int, 'SEQ': str,
        'QUAL': str
    }
    n_fixed_cols = len(fixed_cols)
    cols = fixed_cols + [
        'OPT{}'.format(i)
        for i in range(max(0, (line.count('\t') + 1) - n_fixed_cols))
    ]
    col_dtypes = {k: (fixed_col_dtypes.get(k) or str) for k in cols}
    return pd.read_csv(
        io.StringIO(line), header=None, names=cols, dtype=col_dtypes, sep='\t'
    ).to_dict()
