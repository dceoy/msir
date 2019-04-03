#!/usr/bin/env python

import bz2
from collections import OrderedDict
import gzip
import io
from itertools import chain, product
import logging
import os
import subprocess
from Bio import SeqIO
import pandas as pd
from ..df.beddf import BedDataFrame
from .helper import fetch_abspath, fetch_executable, print_log, \
    run_and_parse_subprocess


def read_fasta(path):
    p = fetch_abspath(path=path)
    if p.endswith('.gz'):
        f = gzip.open(p, 'rt')
    elif p.endswith('.bz2'):
        f = bz2.open(p, 'rt')
    else:
        f = open(p, 'r')
    records = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
    f.close()
    return records


def read_bed(path):
    return BedDataFrame(path=fetch_abspath(path=path)).load_and_output_df()


def iterate_unique_repeat_units(max_unit_len=6, bases='ACGT'):
    patterns_by_ul = []
    for ul in range(1, max_unit_len + 1):
        if ul == 1:
            patterns_by_ul.append(list(bases))
        else:
            unit_tuples_seen = []
            for i in range(1, ul):
                if ul % i == 0:
                    unit_tuples_seen.extend([
                        tuple(s * int(ul / i)) for s in patterns_by_ul[i - 1]
                    ])
            patterns_by_ul.append([
                ''.join(t) for t in product(bases, repeat=ul)
                if t not in unit_tuples_seen
            ])
    return chain.from_iterable(patterns_by_ul)


def convert_bed_line_to_sam_region(bedline):
    return '{0}:{1}-{2}'.format(
        bedline['chrom'], bedline['chromStart'] + 1,  bedline['chromEnd'] + 1
    )


def validate_or_prepare_bam_indexes(bam_paths, index_bam=False, n_proc=8,
                                    samtools_path=None):
    logger = logging.getLogger(__name__)
    invalid_exts = [p for p in bam_paths if not p.endswith(('.bam', '.cram'))]
    if invalid_exts:
        raise RuntimeError('invalid BAM/CRAM paths: {}'.format(invalid_exts))
    else:
        bai_abspaths = {
            p: fetch_abspath(p + '.bai' if p.endswith('.bam') else p + '.crai')
            for p in bam_paths
        }
        bam_paths_without_bai = [
            k for k, v in bai_abspaths.items() if not os.path.isfile(v)
        ]
        if not bam_paths_without_bai:
            pass
        elif index_bam:
            print_log('Index BAM/CRAM files.')
            samtools = samtools_path or fetch_executable('samtools')
            for p in bam_paths_without_bai:
                args = [samtools, 'index', '-@', str(n_proc), fetch_abspath(p)]
                logger.debug('args: {}'.format(args))
                subprocess.run(args=args, check=True)
                print(p + '.bai', flush=True)
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
        (
            samtools, 'view', *options, fetch_abspath(bam_path),
            '{0}:{1}-{1}'.format(rname, p)
        ) for p in [start_pos, end_pos]
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
    ).iloc[0].to_dict(into=OrderedDict)
