#!/usr/bin/env python

import logging
import os
import subprocess


def validate_files_and_dirs(files=[], dirs=[]):
    if files:
        file_nf = [f for f in files if not os.path.isfile(fetch_abspath(f))]
        if file_nf:
            raise FileNotFoundError(
                'File not found: {}'.format(', '.join(file_nf))
            )
    if dirs:
        dir_nf = [d for d in dirs if not os.path.isdir(fetch_abspath(d))]
        if dir_nf:
            raise FileNotFoundError(
                'Directory not found: {}'.format(', '.join(dir_nf))
            )


def fetch_abspath(path):
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def fetch_executable(cmd):
    executables = [
        cp for cp in [
            os.path.join(p, cmd) for p in str.split(os.environ['PATH'], ':')
        ] if os.access(cp, os.X_OK)
    ]
    return executables[0] if executables else None


def print_log(message, prompt='>>>', end=os.linesep):
    print('{0}\t{1}'.format(prompt, message), end=end, flush=True)


def run_and_parse_subprocess(args, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, **kwargs):
    logger = logging.getLogger(__name__)
    logger.debug('args: {}'.format(args))
    with subprocess.Popen(args=args, stdout=stdout, stderr=stderr,
                          **kwargs) as p:
        for line in p.stdout:
            yield line.decode('utf-8')
        outs, errs = p.communicate()
        if p.returncode != 0:
            logger.error(
                'STDERR from subprocess `{0}`:{1}{2}'.format(
                    p.args, os.linesep, errs.decode('utf-8')
                )
            )
            raise subprocess.CalledProcessError(
                returncode=p.returncode, cmd=p.args, output=outs, stderr=errs
            )
