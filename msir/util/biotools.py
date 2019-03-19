#!/usr/bin/env python

import bz2
import gzip
import io
import logging
import os
import subprocess
from Bio import SeqIO
import pandas as pd
from ..util.helper import fetch_executable, print_log, run_and_parse_subprocess


def convert_bed_line_to_sam_region(bedline):
    return '{0}:{1}-{2}'.format(
        bedline['chrom'], bedline['chromStart'] + 1,  bedline['chromEnd'] + 1
    )


def read_fasta(fa_path):
    if fa_path.endswith('.gz'):
        f = gzip.open(fa_path, 'rt')
    elif fa_path.endswith('.bz2'):
        f = bz2.open(fa_path, 'rt')
    else:
        f = open(fa_path, 'r')
    records = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
    f.close()
    return records


def validate_or_prepare_bam_indexes(bam_paths, index_bam=False, n_proc=8,
                                    samtools_path=None):
    logger = logging.getLogger(__name__)
    invalid_exts = [p for p in bam_paths if not p.endswith(('.bam', '.cram'))]
    if invalid_exts:
        raise RuntimeError('invalid BAM/CRAM paths: {}'.format(invalid_exts))
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
                args = [samtools, 'index', '-@', str(n_proc), p]
                logger.debug('args: {}'.format(args))
                subprocess.run(args=args, check=True)
        else:
            raise FileNotFoundError(
                'BAM/CRAM indexes not found: {}'.format(
                    list(bam_paths_without_bai)
                )
            )


def view_bam_lines_including_region(bam_path, rname, start_pos, end_pos,
                                    options=[], samtools_path=None,
                                    parse_lines=True):
    samtools = samtools_path or fetch_executable('samtools')
    args_list = [
        (samtools, 'view', *options, bam_path, '{0}:{1}-{1}'.format(rname, p))
        for p in [start_pos, end_pos]
    ]
    end_pos_set = set(run_and_parse_subprocess(args=args_list[1]))
    for r in run_and_parse_subprocess(args=args_list[0]):
        if r not in end_pos_set:
            pass
        elif parse_lines:
            yield _parse_sam_line(line=r)
        else:
            yield r


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
    ).iloc[0].to_dict()
