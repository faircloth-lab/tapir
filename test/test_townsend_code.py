#!/usr/bin/env python
# encoding: utf-8

"""
test_townsend_code.py

Created by Brant Faircloth on 30 August 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import numpy
import unittest
from models_and_rates import *

import pdb

class TestTransform(unittest.TestCase):
    def setUp(self):
        self.expected = numpy.load('test/test-data/test-parsed-rates.npy')

    def test_uncorrected_site_rates(self):
        observed = parse_site_rates('test/test-data/test-uniform-draw-weights.rates.json')
        assert observed.all() == self.expected.all()

    def test_rate_correction(self):
        observed = parse_site_rates('test/test-data/test-uniform-draw-weights.rates.json', 10.)
        corrected_expected = self.expected/10.
        assert observed.all() == corrected_expected.all()

    def tearDown(self):
        pass

class TestTownsendAndIntegrationAgainstR(unittest.TestCase):
    def setUp(self):
        self.rates = parse_site_rates('test/test-data/test-uniform-draw-weights.rates.json')
        self.time = numpy.array(range(0, 100 + 1, 1))

    def test_townsend_computation_against_R(self):
        townsend = get_townsend_pi(self.time, self.rates)
        assert townsend.all() == numpy.load('test/test-data/test-R-townsend-output.npy').all()

    def test_integration_against_R(self):
        vec_integrate = vectorize(integrate)
        # scipy.integrate returns tuple of (integral, upper-error-bound)
        integral, error = vec_integrate(30, 50, self.rates)
        assert integral.all() == numpy.load('test/test-data/test-30-50-integral.npy').all()

    def tearDown(self):
        pass

class TestTreeAdjustment(unittest.TestCase):
    def setUp(self):
        self.expected_tree = \
            dendropy.Tree.get_from_string('(danRer6:1.74,(oryLat2:1,(gasAcu1:0.93,(fr2:0.37,tetNig2:0.37):0.56):0.07):0.74)', schema='newick')

    def runTest(self):
        depth, correction, self.observed = correct_branch_lengths('test/test-data/Euteleost.tree', 'newick')
        observed_tree = dendropy.Tree.get_from_path(self.observed, 'newick')
        #pdb.set_trace()
        assert observed_tree.as_string('newick') == \
            self.expected_tree.as_string('newick')

    def tearDown(self):
        os.remove(self.observed)

if __name__ == '__main__':
    unittest.main()