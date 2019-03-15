#!/usr/bin/env python
#
# Pandas-based Data Frame Handlers DNA-sequencing
# https://github.com/dceoy/pandna

from abc import ABCMeta, abstractmethod
import os
import subprocess
import pandas as pd


class BaseBioDataFrame(object, metaclass=ABCMeta):
    def __init__(self, path, supported_exts=[]):
        if os.path.isfile(path):
            self.path = path
        else:
            raise BioDataFrameError('file not found: {}'.format(path))
        hit_exts = [x for x in supported_exts if path.endswith(x)]
        if supported_exts and not hit_exts:
            raise BioDataFrameError('invalid file extension: {}'.format(path))
        self.df = pd.DataFrame()

    @abstractmethod
    def load(self):
        pass

    def load_and_output_df(self):
        self.load()
        return self.df

    def write_df(self, path, mode='w', **kwargs):
        if self.header:
            with open(path, mode=mode) as f:
                for h in self.header:
                    f.write(h + os.linesep)
        self.df.to_csv(path, mode=('a' if self.header else 'w'), **kwargs)

    @staticmethod
    def run_and_parse_subprocess(args, stdout=subprocess.PIPE, **kwargs):
        with subprocess.Popen(args=args, stdout=stdout, **kwargs) as p:
            for line in p.stdout:
                yield line.decode('utf-8')
            outs, errs = p.communicate()
            if p.returncode == 0:
                pass
            else:
                raise subprocess.CalledProcessError(
                    returncode=p.returncode, cmd=p.args, output=outs,
                    stderr=errs
                )


class BioDataFrameError(RuntimeError):
    pass
