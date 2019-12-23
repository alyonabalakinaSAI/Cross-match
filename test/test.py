#!/usr/bin/python3

import os
import unittest
from astropy.io import ascii
from astropy.utils.diff import report_diff_values

from CrossMatch import crossmatch


class CrossmatchTestCase(unittest.TestCase):
    def test_match(self):
        testfile = os.path.dirname(os.path.abspath(__file__))
        crossmatch.main('-d ' + testfile + '/crossmatch/testcatalog')
        cat1 = ascii.read('cat1.csv', delimiter=',', guess=False, fast_reader=False)
        cat2 = ascii.read('cat2.csv', delimiter=',', guess=False, fast_reader=False)
        table = crossmatch(cat1, cat2, 'RA1', 'DEC1', 'Z1', 'RA2', 'DEC2', 'Z2', 20, use_z=True)
        final_table = table['NAME1', 'NAME2', 'RA1', 'RA2', 'DEC1', 'DEC2', 'DELTA_Z']
        final_table.write('TEST', format='ascii')
        cat3 = ascii.read('cat3.csv')

        self.assertEqual(True, report_diff_values(cat3, final_table))


if __name__ == '__main__':
    unittest.main()
