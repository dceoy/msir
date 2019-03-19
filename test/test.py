#!/usr/bin/env python

import unittest
from msir.call.identifier import extract_longest_repeat_df, \
    _make_str_regex_dict


class TandemRepeats(unittest.TestCase):
    """Tandem repeat sequence
    """
    reads = {
        'TTTTGCAGAGTACAAGAGTG': {
            'repeat_unit': 'AG', 'repeat_unit_size': 2, 'repeat_seq_size': 4,
            'repeat_times': 2, 'repeat_start': 6, 'repeat_end': 10
        },
        'GTTGGGAAAAAAAAAAATTG': {
            'repeat_unit': 'A', 'repeat_unit_size': 1, 'repeat_seq_size': 11,
            'repeat_times': 11, 'repeat_start': 6, 'repeat_end': 17
        },
        'TATTTTATTATTATTATTAG': {
            'repeat_unit': 'TTA', 'repeat_unit_size': 3, 'repeat_seq_size': 15,
            'repeat_times': 5, 'repeat_start': 4, 'repeat_end': 19
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
