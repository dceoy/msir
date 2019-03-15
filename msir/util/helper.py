#!/usr/bin/env python

from concurrent.futures import as_completed
import logging
import os
from pprint import pformat
import subprocess
import traceback


def validate_files_and_dirs(files=list(), dirs=list()):
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
    """Fetch an executable command path
    Args:
        cmd (str): command name
    """
    executables = [
        cp for cp in [
            os.path.join(p, cmd) for p in str.split(os.environ['PATH'], ':')
        ] if os.access(cp, os.X_OK)
    ]
    return executables[0] if executables else None


def run_and_parse_subprocess(args, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, **kwargs):
    with subprocess.Popen(args=args, stdout=stdout, stderr=stderr,
                          **kwargs) as p:
        for line in p.stdout:
            yield line.decode('utf-8')
        outs, errs = p.communicate()
        if p.returncode == 0:
            pass
        else:
            logger = logging.getLogger(__name__)
            logger.error(
                'STDERR from subprocess `{0}`:{1}{2}'.format(
                    p.args, os.linesep, errs.decode('utf-8')
                )
            )
            raise subprocess.CalledProcessError(
                returncode=p.returncode, cmd=p.args, output=outs, stderr=errs
            )


def wait_futures(futures, ignore_error=False, timeout=604800):
    logger = logging.getLogger(__name__)
    future_results = []
    for f in as_completed(futures):
        try:
            res = f.result(timeout=timeout)
        except Exception as e:
            logger.error(os.linesep + traceback.format_exc())
            if ignore_error:
                logger.debug('Ignore a future error.')
            else:
                canceled_fs = [f for f in futures if f.cancel()]
                if canceled_fs:
                    logger.debug('Canceled futures: {}'.format(canceled_fs))
                raise e
        else:
            logger.info('Complete a future: {0} => {1}'.format(f, res))
            future_results.append(res)
        finally:
            logger.debug('finished future:' + os.linesep + pformat(vars(f)))
    logger.debug('future_results: {}'.format(future_results))
