#!/usr/bin/env python

import unittest
from msir.call.identifier import extract_longest_repeat_df, \
    _make_str_regex_dict


class TandemRepeats(unittest.TestCase):
    """Tandem repeat sequence
    """
    reads = {
        'TTTTGCAGAGTACAAGAGTG': {
            'repeat_unit': 'T', 'repeat_unit_size': 1, 'repeat_times': 4,
            'repeat_start': 0, 'repeat_end': 3
        },
        'GTTGGGAAAAAAAAAAATTG': {
            'repeat_unit': 'A', 'repeat_unit_size': 1, 'repeat_times': 11,
            'repeat_start': 6, 'repeat_end': 16
        },
        'TATTTTATTATTATTATTAG': {
            'repeat_unit': 'TTA', 'repeat_unit_size': 3, 'repeat_times': 5,
            'repeat_start': 4, 'repeat_end': 18
        }
    }

    def test_extract_longest_repeat_df(self):
        """tandem repeat count from sequences
        """
        red = _make_str_regex_dict(max_unit_len=3)
        for s, d in self.reads.items():
            df = extract_longest_repeat_df(sequence=s, regex_dict=red)
            for k, v in d.items():
                self.assertEqual(v, df[k].iloc[0])


if __name__ == '__main__':
    unittest.main()
