#!/usr/bin/env python
"""
Tandem repeat analyzer for microsatellite instability detection by DNA-seq

Usage:
    msir extract [--debug] [--bed=<path>] [--out-dir=<path>] [--index-bam]
                 [--samtools=<path>] [--tsv] [--processes=<int>] <bam>...
    msir -h|--help
    msir -v|--version

Options:
    -h, --help          Print help and exit
    -v, --version       Print version and exit
    --debug             Execute a command with debug messages
    --bed=<path>        Pass an input BED file [default: msir.bed]
    --out-dir=<path>    Pass an output directory [default: .]
    --index-bam         Index BAM or CRAM files if required
    --samtools=<path>   Pass a path to samtools command
    --tsv               Write results into TSV files (CSV is the default)
    --processes=<int>   Limit max cores for multiprocessing

Arguments:
    <bam>               Path to an input BAM/CRAM file

Commands:
    extract             Extract tandem repeats in reads
"""

import logging
import os
from docopt import docopt
from .. import __version__
from ..call.tandemrepeat import extract_tandem_repeats_in_bams


def main():
    args = docopt(__doc__, version='msir {}'.format(__version__))
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=(
            logging.DEBUG if args['--debug']
            else (logging.INFO if args['--info'] else logging.WARNING)
        )
    )
    logging.debug('args:{0}{1}'.format(os.linesep, args))
    if args['extract']:
        extract_tandem_repeats_in_bams(
            bam_paths=args['<bam>'], bed_path=args['--bed'],
            out_dir_path=args['--out-dir'], output_tsv=args['--tsv'],
            index_bam=args['--index-bam'], samtools=args['--samtools'],
            processes=args['--processes']
        )
