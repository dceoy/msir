#!/usr/bin/env python
#
# Pandas-based Data Frame Handlers DNA-sequencing
# https://github.com/dceoy/pandna

import io
import logging
import re
import pandas as pd
from .basebiodf import BaseBioDataFrame


class SamDataFrame(BaseBioDataFrame):
    def __init__(self, path, samtools='samtools', n_thread=1):
        super().__init__(path=path, supported_exts=['.sam', '.bam', '.cram'])
        self.__logger = logging.getLogger(__name__)
        self.__samtools = samtools
        self.__n_th = n_thread
        self.__fixed_cols = [
            'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT',
            'TLEN', 'SEQ', 'QUAL'
        ]
        self.__fixed_col_dtypes = {
            'QNAME': str, 'FLAG': int, 'RNAME': str, 'POS': int, 'MAPQ': int,
            'CIGAR': str, 'RNEXT': str, 'PNEXT': int, 'TLEN': int, 'SEQ': str,
            'QUAL': str
        }
        self.__detected_cols = []
        self.__detected_col_dtypes = {}
        self.header = []

    def load(self):
        if self.path.endswith('.sam'):
            with open(self.path, 'r') as f:
                for s in f:
                    self._load_sam_line(string=s)
        else:
            th_args = (['-@', str(self.__n_th)] if self.__n_th > 1 else [])
            args = [self.__samtools, 'view', *th_args, '-h', self.path]
            for s in self.run_and_parse_subprocess(args=args):
                self._load_sam_line(string=s)
        self.df = self.df.reset_index(drop=True)

    def _load_sam_line(self, string):
        if re.match(r'@[A-Z]{1}', string):
            self.header.append(string.strip())
        else:
            if not self.__detected_cols:
                n_fixed_cols = len(self.__fixed_cols)
                n_detected_cols = string.count('\t') + 1
                self.__detected_cols = self.__fixed_cols + (
                    [
                        'OPT{}'.format(i)
                        for i in range(n_detected_cols - n_fixed_cols)
                    ] if n_detected_cols > n_fixed_cols else []
                )
                self.__detected_col_dtypes = {
                    k: (self.__fixed_col_dtypes.get(k) or str)
                    for k in self.__detected_cols
                }
            self.df = self.df.append(
                pd.read_table(
                    io.StringIO(string), header=None,
                    names=self.__detected_cols,
                    dtype=self.__detected_col_dtypes
                )
            )


class SamtoolsFlagstatDataFrame(BaseBioDataFrame):
    def __init__(self, path):
        super().__init__(path=path)
        self.__cols = ['qc_passed', 'qc_failed', 'read']
        self.__col_dtypes = {'read': str, 'qc_passed': int, 'qc_failed': int}

    def load(self):
        with open(self.path, 'r') as f:
            for s in f:
                self._load_samtools_flagstat_line(string=s)
        self.df = self.df.reset_index(drop=True)

    def _load_samtools_flagstat_line(self, string):
        self.df = self.df.append(
            pd.read_table(
                io.StringIO(
                    string.replace(' + ', '\t', 1).replace(' ', '\t', 1)
                ),
                header=None, names=self.__cols, dtype=self.__col_dtypes
            )[['read', 'qc_passed', 'qc_failed']]
        )
