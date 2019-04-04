#!/usr/bin/env python

from collections import OrderedDict
from concurrent.futures import as_completed, ProcessPoolExecutor
import logging
import os
import pandas as pd
from ..util.helper import fetch_abspath, print_log, validate_files_and_dirs
from ..util.biotools import convert_bed_line_to_sam_region, \
    validate_or_prepare_bam_indexes, view_bam_lines_including_region
from .identifier import compile_str_regex, extract_longest_repeat_df


def detect_tandem_repeats_in_reads(bam_paths, trunit_tsv_path, obs_tsv_path,
                                   index_bam=False, append_read_seq=False,
                                   cut_end_len=10, samtools=None, n_proc=8):
    validate_files_and_dirs(files=[trunit_tsv_path, *bam_paths])
    validate_or_prepare_bam_indexes(
        bam_paths=bam_paths, index_bam=index_bam, n_proc=n_proc,
        samtools_path=samtools
    )
    print_log('Load repeat units data:\t{}'.format(trunit_tsv_path))
    df_ru = pd.read_csv(trunit_tsv_path, sep='\t')
    print_log('Compile regular expression patterns:', end='')
    regex_patterns = _compile_repeat_unit_regex_patterns_from_df(df=df_ru)
    print('\t{}'.format(len(regex_patterns)), flush=True)
    tsv_abspath = fetch_abspath(obs_tsv_path)
    for i, p in enumerate(bam_paths):
        print_log('Detect tandem repeats within reads:\t{}'.format(p))
        df_obs = _extract_repeats_within_reads(
            bam_path=p, regex_patterns=regex_patterns, df_ru=df_ru,
            cut_end_len=int(cut_end_len), append_read_seq=append_read_seq,
            samtools=samtools, n_proc=n_proc
        )
        print_log('Write repeat counts data:\t{}'.format(obs_tsv_path))
        df_obs.to_csv(
            tsv_abspath, header=(not i), mode=('a' if i else 'w'),
            sep=(',' if tsv_abspath.endswith('.csv') else '\t')
        )
    print_log('All the processes done.')


def _compile_repeat_unit_regex_patterns_from_df(df):
    return {
        s: compile_str_regex(repeat_unit=s, min_rep_times=1)
        for s in set(df['repeat_unit'])
    }


def _extract_repeats_within_reads(bam_path, regex_patterns, df_ru,
                                  append_read_seq, cut_end_len, samtools,
                                  n_proc=8):
    logger = logging.getLogger(__name__)
    ppx = ProcessPoolExecutor(max_workers=n_proc)
    fs = [
        ppx.submit(
            _count_repeats_in_reads, bam_path, line, id, regex_patterns,
            append_read_seq, cut_end_len, samtools
        ) for id, line in df_ru.iterrows()
    ]
    try:
        df_obs = pd.concat(
            [f.result() for f in as_completed(fs)], sort=False
        ).sort_values(by='id').drop(columns='id')
    except Exception as e:
        [f.cancel() for f in fs]
        ppx.shutdown(wait=False)
        raise e
    else:
        ppx.shutdown(wait=True)
    finally:
        logger.debug('df_obs:{0}{1}'.format(os.linesep, df_obs))
    return df_obs


def _count_repeats_in_reads(bam_path, tsvline, id, regex_patterns,
                            append_read_seq, cut_end_len, samtools):
    logger = logging.getLogger(__name__)
    bed_cols = ['chrom', 'chromStart', 'chromEnd']
    sam_cols = [
        'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT',
        'TLEN', *(['SEQ', 'QUAL'] if append_read_seq else [])
    ]
    ru = tsvline['repeat_unit']
    region = convert_bed_line_to_sam_region(tsvline)
    logger.debug('region: {}'.format(region))
    view_args = {
        'bam_path': fetch_abspath(bam_path), 'rname': tsvline['chrom'],
        'start_pos': (tsvline['repeat_start'] + 1 - cut_end_len),
        'end_pos': (tsvline['repeat_end'] + cut_end_len),
        'samtools_path': samtools
    }
    line_df_list = []
    for d in view_bam_lines_including_region(**view_args):
        df = extract_longest_repeat_df(
            sequence=d['SEQ'], regex_patterns={ru: regex_patterns[ru]},
            cut_end_len=cut_end_len, start_pos=d['POS']
        )
        if df.size:
            sam_od = OrderedDict([(k, d[k]) for k in sam_cols])
            line_df_list.append(df.drop(columns='repeat_seq').assign(**sam_od))
    if line_df_list:
        df_obs = pd.concat(line_df_list, sort=False).rename(
            columns={
                'repeat_start': 'observed_repeat_start',
                'repeat_end': 'observed_repeat_end',
                'repeat_seq_length': 'observed_repeat_seq_length',
                'repeat_times': 'observed_repeat_times'
            }
        ).assign(
            id=id, data_path=bam_path, sam_region=region,
            ref_repeat_times=tsvline['repeat_times'],
            **{k: tsvline[k] for k in bed_cols}
        ).set_index([
            'data_path', *bed_cols, 'sam_region', 'repeat_unit',
            'repeat_unit_length', 'ref_repeat_times'
        ]).sort_index()
        logger.debug('df_obs:{0}{1}'.format(os.linesep, df_obs))
        _print_state_line(region=region, df_obs=df_obs, bam_path=bam_path)
    else:
        df_obs = pd.DataFrame()
    return df_obs


def _print_state_line(region, df_obs, bam_path):
    bam_name = os.path.basename(bam_path)
    ru = df_obs.reset_index()['repeat_unit'].iloc[0]
    df_hist = df_obs['observed_repeat_times'].value_counts().to_frame(name='c')
    for i, r in df_hist.iterrows():
        line = '  {0}\t{1:<25}\t{2:<10}\t{3}'.format(
            bam_name, region, '{0}x{1}'.format(i, ru), r['c']
        )
        print(line, flush=True)
