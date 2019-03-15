#!/usr/bin/env python

import logging
from multiprocessing import cpu_count
from ..df.beddf import BedDataFrame
from ..util.helper import validate_files_and_dirs, fetch_executable


def extract_tandem_repeats_in_sams(bed_path, sam_paths, out_dir_path,
                                   output_tsv=False, samtools=None,
                                   processes=None):
    validate_files_and_dirs(files=[bed_path, *sam_paths], dirs=[out_dir_path])
    common_args = {
        'bed_path': bed_path, 'out_dir_path': out_dir_path,
        'output_tsv': output_tsv,
        'samtools': (samtools or fetch_executable('samtools')),
        'n_proc': int(processes or cpu_count())
    }
    [_extract_tandem_repeats(sam_path=p, **common_args) for p in sam_paths]


def _extract_tandem_repeats(sam_path, bed_path, out_dir_path, output_tsv,
                            samtools, n_proc):
    logger = logging.getLogger(__name__)
    print('msir>\tExtract repeats in:\t{}'.format(sam_path), flush=True)
