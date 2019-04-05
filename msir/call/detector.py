#!/usr/bin/env python

from collections import OrderedDict
from concurrent.futures import as_completed, ProcessPoolExecutor
import logging
import os
import traceback
import pandas as pd
from ..util.helper import fetch_abspath, print_log, validate_files_and_dirs
from ..util.biotools import convert_bed_line_to_sam_region, \
    validate_or_prepare_bam_indexes, view_bam_lines_including_region
from .identifier import compile_str_regex, extract_longest_repeat_df


def detect_tandem_repeats_in_reads(bam_paths, trunit_tsv_path, obs_tsv_path,
                                   index_bam=False, append_read_seq=False,
                                   samtools=None, n_proc=8):
    validate_files_and_dirs(files=[trunit_tsv_path, *bam_paths])
    validate_or_prepare_bam_indexes(
        bam_paths=bam_paths, index_bam=index_bam, n_proc=n_proc,
        samtools_path=samtools
    )
    print_log('Load repeat units data:\t{}'.format(trunit_tsv_path))
    df_ru = pd.read_csv(trunit_tsv_path, sep='\t')
    print_log('Compile regular expression patterns:', end='')
    regex_dict = _compile_repeat_unit_regex_patterns_from_df(df=df_ru)
    print('\t{}'.format(len(regex_dict)), flush=True)
    tsv_abspath = fetch_abspath(obs_tsv_path)
    if os.path.exists(tsv_abspath):
        os.remove(tsv_abspath)
    for i, p in enumerate(bam_paths):
        print_log('Detect tandem repeats within reads:\t{}'.format(p))
        df_obs = _extract_repeats_within_reads(
            bam_path=p, regex_dict=regex_dict, df_ru=df_ru,
            append_read_seq=append_read_seq, samtools=samtools, n_proc=n_proc
        )
        if df_obs.size:
            print_log('Write observed repeat data:\t{}'.format(obs_tsv_path))
            df_obs.to_csv(
                tsv_abspath, header=(not os.path.exists(tsv_abspath)),
                mode='a', sep=(',' if tsv_abspath.endswith('.csv') else '\t')
            )
        else:
            print_log('No repeats were detected:\t{}'.format(p))
    print_log('All the processes done.')


def _compile_repeat_unit_regex_patterns_from_df(df):
    return {
        i: {
            'patterns': {
                r['repeat_unit']: compile_str_regex(
                    repeat_unit=r['repeat_unit'], min_rep_times=1,
                    left_seq=r['left_seq'], right_seq=r['right_seq']
                )
            },
            'left_seq': r['left_seq'], 'right_seq': r['right_seq']
        } for i, r in df.iterrows()
    }


def _extract_repeats_within_reads(bam_path, regex_dict, df_ru, append_read_seq,
                                  samtools, n_proc=8):
    logger = logging.getLogger(__name__)
    ppx = ProcessPoolExecutor(max_workers=n_proc)
    fs = [
        ppx.submit(
            _extract_repeats_at_region, bam_path, line, id, regex_dict,
            append_read_seq, samtools
        ) for id, line in df_ru.iterrows()
    ]
    try:
        df_obs = pd.concat(
            [f.result() for f in as_completed(fs)], sort=False
        )
    except Exception as e:
        logger.error(os.linesep + traceback.format_exc())
        ppx.shutdown(wait=False)
        raise e
    else:
        logger.debug('df_obs:{0}{1}'.format(os.linesep, df_obs))
        ppx.shutdown(wait=True)
    if df_obs.size:
        return df_obs.sort_values(by='bed_id').drop(columns='bed_id')
    else:
        return df_obs


def _extract_repeats_at_region(bam_path, tsvline, bed_id, regex_dict,
                               append_read_seq, samtools):
    logger = logging.getLogger(__name__)
    bed_cols = ['chrom', 'chromStart', 'chromEnd']
    sam_cols = [
        'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT',
        'TLEN', *(['SEQ', 'QUAL'] if append_read_seq else [])
    ]
    region = convert_bed_line_to_sam_region(tsvline)
    logger.debug('region: {}'.format(region))
    view_args = {
        'bam_path': fetch_abspath(bam_path), 'rname': tsvline['chrom'],
        'start_pos': (tsvline['repeat_start'] + 1 - len(tsvline['left_seq'])),
        'end_pos': (tsvline['repeat_end'] + len(tsvline['right_seq'])),
        'samtools_path': samtools
    }
    region_df_list = []
    regex_patterns = regex_dict[bed_id]
    logger.debug('regex_patterns:\t{}'.format(regex_patterns))
    region_df_list_raw = [
        extract_longest_repeat_df(
            sequence=b['SEQ'], regex_patterns=regex_patterns, min_rep_len=1,
            flanking_len=0, start_pos=b['POS']
        ).pipe(
            lambda d: d.drop(
                columns=['start_x', 'end_x', 'repeat_seq']
            ).assign(
                **OrderedDict([(k, b[k]) for k in sam_cols])
            ) if d.size else d
        ) for b in view_bam_lines_including_region(**view_args)
    ]
    region_df_list = [d for d in region_df_list_raw if d.size]
    if not region_df_list:
        return pd.DataFrame()
    else:
        df_region = pd.concat(region_df_list, sort=False).rename(
            columns={
                'repeat_start': 'observed_repeat_start',
                'repeat_end': 'observed_repeat_end',
                'repeat_seq_length': 'observed_repeat_seq_length',
                'repeat_times': 'observed_repeat_times'
            }
        ).assign(
            bed_id=bed_id, sam_path=bam_path, sam_region=region,
            referenced_repeat_times=tsvline['repeat_times'],
            **{k: tsvline[k] for k in bed_cols}
        ).set_index([
            'sam_path', *bed_cols, 'sam_region', 'repeat_unit',
            'repeat_unit_length', 'referenced_repeat_times', 'left_seq',
            'right_seq'
        ]).sort_index()
        logger.debug('df_region:{0}{1}'.format(os.linesep, df_region))
        _print_state_line(region=region, df=df_region, bam_path=bam_path)
        return df_region


def _print_state_line(region, df, bam_path):
    bam_name = os.path.basename(bam_path)
    ru = df.reset_index()['repeat_unit'].iloc[0]
    df_hist = df['observed_repeat_times'].value_counts().to_frame(name='c')
    for i, r in df_hist.iterrows():
        line = '  {0}\t{1:<25}\t{2:<10}\t{3}'.format(
            bam_name, region, '{0}x{1}'.format(i, ru), r['c']
        )
        print(line, flush=True)
