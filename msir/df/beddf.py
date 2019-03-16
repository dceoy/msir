#!/usr/bin/env python
#
# Pandas-based Data Frame Handlers DNA-sequencing
# https://github.com/dceoy/pandna

import io
import pandas as pd
from .basebiodf import BaseBioDataFrame


class BedDataFrame(BaseBioDataFrame):
    def __init__(self, path, opt_cols=[]):
        super().__init__(path=path, supported_exts=['.bed', '.txt', '.tsv'])
        self.__fixed_cols = ['chrom', 'chromStart', 'chromEnd']
        self.__opt_cols = opt_cols or [
            'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb',
            'blockCount', 'blockSizes', 'blockStarts'
        ]
        self.__fixed_col_dtypes = {
            'chrom': str, 'chromStart': int, 'chromEnd': int, 'name': str,
            'score': int, 'strand': str, 'thickStart': int, 'thickEnd': int,
            'itemRgb': str, 'blockCount': int, 'blockSizes': int,
            'blockStarts': int
        }
        self.__detected_cols = []
        self.__detected_col_dtypes = {}
        self.header = []

    def load(self):
        with open(self.path, 'r') as f:
            for s in f:
                self._load_bed_line(string=s)
        self.df = self.df.reset_index(drop=True)

    def _load_bed_line(self, string):
        if string.startswith(('browser', 'track')):
            self.header.append(string.strip())
        else:
            if not self.__detected_cols:
                self.__detected_cols = [
                    *self.__fixed_cols, *self.__opt_cols
                ][:(string.count('\t') + 1)]
                self.__detected_col_dtypes = {
                    k: (self.__fixed_col_dtypes.get(k) or str)
                    for k in self.__detected_cols
                }
            self.df = self.df.append(
                pd.read_csv(
                    io.StringIO(string), header=None,
                    names=self.__detected_cols,
                    dtype=self.__detected_col_dtypes, sep='\t'
                )
            )
