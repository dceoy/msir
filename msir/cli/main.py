#!/usr/bin/env python
"""
Tandem repeat analyzer for microsatellite instability detection by DNA-seq

Usage:
    msir extract [--debug] [--bed=<path>] [--out-dir=<path>]
                 [--samtools=<path>] [--tsv] [--processes=<int>] <sam>...
    msir -h|--help
    msir -v|--version

Options:
    -h, --help          Print help and exit
    -v, --version       Print version and exit
    --debug             Execute a command with debug messages
    --bed=<path>        Pass an input BED file [default: msir.bed]
    --out-dir=<path>    Pass an output directory [default: .]
    --samtools=<path>   Pass a path to samtools command
    --tsv               Write results into TSV files (CSV is the default)
    --processes=<int>   Limit max cores for multiprocessing

Arguments:
    <sam>               Path to an input SAM/BAM/CRAM file

Commands:
    extract             Extract tandem repeats in reads
"""

import logging
from docopt import docopt
from .. import __version__
from ..call.tandemrepeat import extract_tandem_repeats_in_sams


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
    logging.debug('args:\n{}'.format(args))
    if args['extract']:
        extract_tandem_repeats_in_sams(
            sam_paths=args['<sam>'], bed_path=args['--bed'],
            out_dir_path=args['--out-dir'], output_tsv=args['--tsv'],
            samtools=args['--samtools'], processes=args['--processes']
        )
