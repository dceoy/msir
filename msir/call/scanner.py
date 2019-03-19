#!/usr/bin/env python

from concurrent.futures import as_completed, ProcessPoolExecutor
import logging
import os
import pandas as pd
from ..util.helper import fetch_abspath, print_log, validate_files_and_dirs
from ..util.biotools import convert_bed_line_to_sam_region, \
    validate_or_prepare_bam_indexes, view_bam_lines_including_region
from .identifier import compile_str_regex, extract_longest_repeat_df


def scan_tandem_repeats_in_reads(bam_paths, ru_tsv_path, result_tsv_path,
                                 index_bam=False, samtools=None,
                                 cut_end_len=10, n_proc=8):
    validate_files_and_dirs(files=[ru_tsv_path, *bam_paths])
    validate_or_prepare_bam_indexes(
        bam_paths=bam_paths, index_bam=index_bam, n_proc=n_proc,
        samtools_path=samtools
    )
    print_log('Load repeat units data:\t{}'.format(ru_tsv_path))
    df_ru = pd.read_csv(ru_tsv_path, sep='\t')
    print_log('Compile regular expression patterns:')
    regex_patterns = _compile_repeat_unit_regex_patterns_from_df(df=df_ru)
    for p in bam_paths:
        print_log('Extract tandem repeats within reads:\t{}'.format(p))
        df_result = _extract_repeats_within_reads(
            bam_path=p, regex_patterns=regex_patterns, df_ru=df_ru,
            cut_end_len=cut_end_len, samtools=samtools, n_proc=n_proc
        )
        print_log('Write results into:\t{}'.format(result_tsv_path))
        _write_result_tsv(df=df_result, tsv_path=result_tsv_path, bam_path=p)
    print_log('All the processes done.')


def _compile_repeat_unit_regex_patterns_from_df(df):
    patterns = {
        s: compile_str_regex(repeat_unit=s, min_rep_times=1)
        for s in set(df['repeat_unit'])
    }
    print('  Patterns: {}'.format(len(patterns)), flush=True)
    return patterns


def _extract_repeats_within_reads(bam_path, regex_patterns, df_ru, cut_end_len,
                                  samtools, n_proc=8):
    ppx = ProcessPoolExecutor(max_workers=n_proc)
    fs = [
        ppx.submit(
            _count_repeats_in_reads, bam_path, line, id, regex_patterns,
            cut_end_len, samtools
        ) for id, line in df_ru.iterrows()
    ]
    df_result = pd.DataFrame()
    try:
        for f in as_completed(fs):
            res = f.result()
            _print_state_line(res=res, bam_path=bam_path)
            df_result = df_result.append(res['df'])
    except Exception as e:
        [f.cancel() for f in fs]
        ppx.shutdown(wait=False)
        raise e
    else:
        ppx.shutdown(wait=True)
        df_result = df_result.set_index('id').sort_index()
    finally:
        print(df_result, flush=True)
    return df_result


def _count_repeats_in_reads(bam_path, tsvline, id, regex_patterns, cut_end_len,
                            samtools):
    logger = logging.getLogger(__name__)
    bed_cols = ['chrom', 'chromStart', 'chromEnd']
    ru = tsvline['repeat_unit']
    region = convert_bed_line_to_sam_region(tsvline)
    logger.debug('region: {}'.format(region))
    view_args = {
        'bam_path': fetch_abspath(bam_path), 'rname': tsvline['chrom'],
        'start_pos': (tsvline['repeat_start'] + 1 - cut_end_len),
        'end_pos': (tsvline['repeat_end'] + cut_end_len),
        'samtools_path': samtools
    }
    df_bam = pd.DataFrame()
    for d in view_bam_lines_including_region(**view_args):
        df = extract_longest_repeat_df(
            sequence=d['SEQ'], regex_patterns={ru: regex_patterns[ru]},
            cut_end_len=cut_end_len, start_pos=d['POS']
        )
        if df.size:
            df_bam = df_bam.append(df)
    logger.debug('df_bam:{0}{1}'.format(os.linesep, df_bam))
    if df_bam.size:
        df_rt = df_bam['repeat_times'].value_counts().to_frame(
            name='observed_repeat_times_count'
        ).reset_index().rename(
            columns={'index': 'observed_repeat_times'}
        ).assign(
            id=id, sam_region=region, ref_repeat_times=tsvline['repeat_times'],
            **{k: tsvline[k] for k in bed_cols}, repeat_unit=ru
        ).sort_values(
            'observed_repeat_times_count', ascending=False
        )[[
            'id', *bed_cols, 'sam_region', 'repeat_unit', 'ref_repeat_times',
            'observed_repeat_times', 'observed_repeat_times_count'
        ]]
    else:
        df_rt = pd.DataFrame()
    return {'region': region, 'df': df_rt}


def _print_state_line(res, bam_path):
    if res['df'].size:
        for i, r in res['df'].iterrows():
            line = '  {0}\t{1:<25}\t{2:<10}\t{3}'.format(
                os.path.basename(bam_path), res['region'],
                '{0}x{1}'.format(r['observed_repeat_times'], r['repeat_unit']),
                r['observed_repeat_times_count']
            )
            print(line, flush=True)


def _write_result_tsv(df, tsv_path, bam_path):
    tsv_abspath = fetch_abspath(tsv_path)
    tsv_exists = os.path.isfile(tsv_abspath)
    df.pipe(
        lambda d: d.assign(data_path=bam_path)[['data_path', *d.columns]]
    ).to_csv(
        tsv_abspath, index=False, header=(not tsv_exists),
        mode=('a' if tsv_exists else 'w'),
        sep=(',' if tsv_path.endswith('.csv') else '\t')
    )
