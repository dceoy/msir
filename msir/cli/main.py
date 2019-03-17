#!/usr/bin/env python
"""
Tandem repeat analyzer for microsatellite instability detection by DNA-seq

Usage:
    msir id [--debug] [--max-unit-len=<int>] [--min-rep-times=<int>]
            [--ex-region-len=<int>] [--processes=<int>] [--unit-tsv=<path>]
            <bed> <fasta>
    msir count [--debug] [--unit-tsv=<path>] [--out-dir=<path>] [--index-bam]
               [--samtools=<path>] [--cut-end-len=<int>] [--csv]
               [--processes=<int>] <bam>...
    msir -h|--help
    msir -v|--version

Options:
    -h, --help              Print help and exit
    -v, --version           Print version and exit
    --debug                 Execute a command with debug messages
    --max-unit-len=<int>    Set a maximum length for repeat units [default: 6]
    --min-rep-times=<int>   Set a minimum repeat times [default: 3]
    --ex-region-len=<int>   Search around extra regions [default: 20]
    --processes=<int>       Limit max cores for multiprocessing
    --unit-tsv=<path>       Set a TSV file of repeat units [default: ru.tsv]
    --out-dir=<path>        Pass an output directory [default: .]
    --index-bam             Index BAM or CRAM files if required
    --samtools=<path>       Set a path to samtools command
    --cut-end-len=<int>     Ignore repeats on ends of reads [default: 10]
    --csv                   Write results with CSV instead of TSV

Arguments:
    <bed>                   Path to a BED file of repetitive regions
    <fasta>                 Path to a reference genome fasta file
    <bam>                   Path to an input BAM/CRAM file

Commands:
    id                      Indentify repeat units from reference sequences
    count                   Extract and count tandem repeats in read sequences
"""

import logging
from multiprocessing import cpu_count
import os
from docopt import docopt
from .. import __version__
from ..call.identifier import identify_repeat_units_on_bed
from ..call.scanner import scan_tandem_repeats_in_reads


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
    logger = logging.getLogger(__name__)
    logger.debug('args:{0}{1}'.format(os.linesep, args))
    n_proc = int(args['--processes'] or cpu_count())
    logger.debug('n_proc: {}'.format(n_proc))
    if args['count']:
        scan_tandem_repeats_in_reads(
            bam_paths=args['<bam>'], ru_tsv_path=args['--unit-tsv'],
            out_dir_path=args['--out-dir'], output_csv=args['--csv'],
            index_bam=args['--index-bam'], samtools=args['--samtools'],
            cut_end_len=args['--cut-end-len'], n_proc=n_proc
        )
    elif args['id']:
        identify_repeat_units_on_bed(
            bed_path=args['<bed>'], genome_fa_path=args['<fasta>'],
            ru_tsv_path=args['--unit-tsv'],
            max_unit_len=int(args['--max-unit-len']),
            min_rep_times=int(args['--min-rep-times']),
            ex_region_len=int(args['--ex-region-len']),
            n_proc=n_proc
        )
