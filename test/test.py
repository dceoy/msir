#!/usr/bin/env python

import unittest
from msir.call.identifier import compile_str_regex, \
    extract_longest_repeat_df
from msir.util.biotools import iterate_unique_repeat_units


class TandemRepeats(unittest.TestCase):
    """Tandem repeat sequence
    """
    repeat_units = {
        1: {'A', 'C', 'G', 'T'},
        2: {'A', 'C', 'G', 'T', 'AC', 'AG', 'AT', 'CA', 'CG', 'CT', 'GA',
            'GC', 'GT', 'TA', 'TC', 'TG'}
    }
    reads = {
        'TTTTGCAGAGAGTACAAGAGG': {
            'repeat_unit': 'AG', 'repeat_unit_length': 2,
            'repeat_seq_length': 6, 'repeat_times': 3, 'repeat_start': 6,
            'repeat_end': 12
        },
        'GTTGGGAAAAAAAAAAATTG': {
            'repeat_unit': 'A', 'repeat_unit_length': 1,
            'repeat_seq_length': 11, 'repeat_times': 11, 'repeat_start': 6,
            'repeat_end': 17
        },
        'TATTTTATTATTATTATTAG': {
            'repeat_unit': 'TTA', 'repeat_unit_length': 3,
            'repeat_seq_length': 15, 'repeat_times': 5, 'repeat_start': 4,
            'repeat_end': 19
        }
    }

    def test_iterate_unique_repeat_units(self, max_unit_len=2):
        for i in range(1, max_unit_len + 1):
            uset = set(iterate_unique_repeat_units(max_unit_len=i))
            self.assertEqual(uset, self.repeat_units[i])

    def test_extract_longest_repeat_df(self):
        """tandem repeat count from sequences
        """
        for s, d in self.reads.items():
            df = extract_longest_repeat_df(
                sequence=s,
                regex_patterns={
                    d['repeat_unit']: compile_str_regex(d['repeat_unit'])
                }
            )
            for k, v in d.items():
                self.assertEqual(v, df[k].iloc[0])


if __name__ == '__main__':
    unittest.main()
