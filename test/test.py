#!/usr/bin/env python

import unittest
from msir.call.identifier import count_repeat_times


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
        'TATTTATTATTATTATTATT': {
            'repeat_unit': 'TTA', 'repeat_unit_size': 3, 'repeat_times': 5,
            'repeat_start': 3, 'repeat_end': 17
        }
    }

    def test_count_repeat_times(self):
        """tandem repeat count from sequences
        """
        for s, d in self.reads.items():
            r = count_repeat_times(sequence=s, repeat_unit=d['repeat_unit'])
            for k, v in d.items():
                self.assertEqual(v, r.get(k))


if __name__ == '__main__':
    unittest.main()
