"""
File: test_base.py
Author: Brant Faircloth

Created by Brant Faircloth on 23 October 2011 16:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: test methods for tapir.base

"""

import shutil
import unittest
import tempfile
from tapir.base import *
from tapir import get_test_files

#import pdb

class TestTransform(unittest.TestCase):

    def setUp(self):
        self.loc = get_test_files()

    def test_is_dir(self):
        # create temp dir
        name = tempfile.mkdtemp()
        assert is_dir(name)
        shutil.rmtree(name)

    def test_get_output_type(self):
        assert get_output_type('test.jpg') == 'jpg'
        # wrong type
        self.assertRaises(AssertionError, get_output_type, 'test.bob')
        # no type
        self.assertRaises(AssertionError, get_output_type, 'test')

    def test_get_list_from_ints(self):
        assert get_list_from_ints('1,2,3') == [1,2,3]

    def test_strings_from_items(self):
        assert get_strings_from_items('1,2,3') == ['1','2','3']

    def test_get_list_from_ranges(self):
        assert get_list_from_ranges('1-2,2-3,3-4') == \
                [[1,2],[2,3],[3,4]]

    def test_get_files_1(self):
        observed = [os.path.basename(i) for i in
                get_files(self.loc,'*.nex,*.nexus')]
        expected = [
                'chr1_918.nex',
                'informativeness_cutoff.nex',
                'test-extension.nexus'
            ]
        assert observed == expected

    def test_get_files_2(self):
        self.assertRaises(IOError, get_files, 'test-data','*.rrwrr')

if __name__ == '__main__':
    unittest.main()
