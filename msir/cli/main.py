#!/usr/bin/env python
"""
Tandem repeat analyzer for microsatellite instability detection by DNA-seq

Usage:
    msir id [--debug] [--unit-tsv=<path>] [--max-unit-len=<int>]
            [--min-rep-times=<int>] [--min-rep-len=<int>]
            [--ex-region-len=<int>] [--processes=<int>] <bed> <fasta>
    msir detect [--debug] [--unit-tsv=<path>] [--obs-tsv=<path>][--index-bam]
                [--samtools=<path>] [--cut-end-len=<int>] [--processes=<int>]
                <bam>...
    msir pipeline [--debug] [--max-unit-len=<int>] [--min-rep-times=<int>]
                  [--min-rep-len=<int>] [--ex-region-len=<int>]
                  [--cut-end-len=<int>] [--index-bam] [--samtools=<path>]
                  [--processes=<int>] [--unit-tsv=<path>] [--obs-tsv=<path>]
                  <bed> <fasta> <bam>...
    msir -h|--help
    msir -v|--version

Options:
    -h, --help              Print help and exit
    -v, --version           Print version and exit
    --debug                 Execute a command with debug messages
    --max-unit-len=<int>    Set a maximum length for repeat units [default: 8]
    --min-rep-times=<int>   Set a minimum repeat times [default: 2]
    --min-rep-len=<int>     Set a minimum length for repeats [default: 5]
    --ex-region-len=<int>   Search around extra regions [default: 20]
    --processes=<int>       Limit max cores for multiprocessing
    --unit-tsv=<path>       Set a TSV of repeat units [default: tr_unit.tsv]
    --obs-tsv=<path>        Set a TSV of observed repeats [default: tr_obs.tsv]
    --index-bam             Index BAM or CRAM if required
    --samtools=<path>       Set a path to samtools command
    --cut-end-len=<int>     Ignore repeats on ends of reads [default: 10]

Arguments:
    <bed>                   Path to a BED file of repetitive regions
    <fasta>                 Path to a reference genome FASTA file
    <bam>                   Path to an input BAM/CRAM file

Commands:
    id                      Indentify repeat units from reference sequences
    detect                  Detect tandem repeats within read sequences
    pipeline                Execute both of the above commands
"""

import logging
from multiprocessing import cpu_count
import os
import signal
from docopt import docopt
from .. import __version__
from ..call.identifier import identify_repeat_units_on_bed
from ..call.detector import detect_tandem_repeats_in_reads


def main():
    args = docopt(__doc__, version='msir {}'.format(__version__))
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=(logging.DEBUG if args['--debug'] else logging.WARNING)
    )
    logger = logging.getLogger(__name__)
    logger.debug('args:{0}{1}'.format(os.linesep, args))
    n_proc = int(args['--processes'] or cpu_count())
    logger.debug('n_proc: {}'.format(n_proc))
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    if args['id'] or args['pipeline']:
        identify_repeat_units_on_bed(
            bed_path=args['<bed>'], genome_fa_path=args['<fasta>'],
            trunit_tsv_path=args['--unit-tsv'],
            max_unit_len=int(args['--max-unit-len']),
            min_rep_times=int(args['--min-rep-times']),
            min_rep_len=int(args['--min-rep-len']),
            ex_region_len=int(args['--ex-region-len']),
            n_proc=n_proc
        )
    if args['detect'] or args['pipeline']:
        detect_tandem_repeats_in_reads(
            bam_paths=args['<bam>'], trunit_tsv_path=args['--unit-tsv'],
            obs_tsv_path=args['--obs-tsv'], index_bam=args['--index-bam'],
            samtools=args['--samtools'],
            cut_end_len=int(args['--cut-end-len']), n_proc=n_proc
        )
