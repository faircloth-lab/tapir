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
from tapir.compute import *
from tapir import get_test_files

import pdb

class TestTransform(unittest.TestCase):
    def setUp(self):
        self.loc = get_test_files()
        self.expected = numpy.load(os.path.join(self.loc,'test-parsed-rates.npy'))

    def test_uncorrected_site_rates(self):
        observed = parse_site_rates(os.path.join(self.loc,'test-uniform-draw-weights.rates.json'),
            test = True)
        assert observed.all() == self.expected.all()

    def test_rate_correction(self):
        observed = parse_site_rates(os.path.join(self.loc,'test-uniform-draw-weights.rates.json'),
                10., test = True)
        corrected_expected = self.expected/10.
        assert observed.all() == corrected_expected.all()

    def tearDown(self):
        pass

class TestTownsendAndIntegrationAgainstPhydesign(unittest.TestCase):
    def setUp(self):
        self.loc = get_test_files()
        self.rates = parse_site_rates(os.path.join(self.loc,'test-uniform-draw-weights.rates.json'), test = True)
        self.time = get_time(0, 174)

    def test_townsend_pi_computation_against_phydesign(self):
        townsend = get_townsend_pi(self.time, numpy.array(self.rates))
        # these are values output from phydesign for the uniform-draw-rates
        # that are part of the json file above, to ensure there is no weird
        # model crap getting in the way
        self.assertAlmostEqual(
            numpy.sum(townsend, axis = 1)[:5].all(),
            numpy.array([0.00000,0.10778,0.14484,0.14616,0.13132,0.11089]).all(),
            3
        )
        expected = [0.03448,0.01111,0.02293]
        for k,v in enumerate([10,20,50]):
            self.assertAlmostEqual(
                numpy.sum(townsend, axis = 1)[v],
                expected[k],
                3
        )

    def test_townsend_pi_integration_against_phydesign(self):
        vec_integrate = vectorize(get_integral_over_times)
        # these are values output from phydesign for the uniform-draw-rates
        # that are part of the json file above, to ensure there is no weird
        # model crap getting in the way
        expected = [0.93453,0.10628,0.05855,0.12638,1.03698,2.08840]
        for k, pair in enumerate(([0,10], [10,15], [15,20], [20,30], [20,70], [20,100])):
            integral, error = vec_integrate(pair[0], pair[1], self.rates)
            #pdb.set_trace()
            self.assertAlmostEqual(sum(integral), expected[k], 4)

    def cleanUp(self):
        pass

class TestTreeAdjustment(unittest.TestCase):
    def setUp(self):
        self.loc = get_test_files()
        self.expected_tree = \
            dendropy.Tree.get_from_string('(danRer6:1.74,(oryLat2:1,(gasAcu1:0.93,(fr2:0.37,tetNig2:0.37):0.56):0.07):0.74)', schema='newick')

    def test_dendropy_tree_adjustment(self):
        depth, correction, self.observed = \
            correct_branch_lengths(os.path.join(self.loc,'Euteleost.tree'), 'newick')
        observed_tree = dendropy.Tree.get_from_path(self.observed, 'newick')
        #pdb.set_trace()
        assert observed_tree.as_string('newick') == \
            self.expected_tree.as_string('newick')

    def tearDown(self):
        os.remove(self.observed)

class TestInformativenessCutoff(unittest.TestCase):
    def setUp(self):
        self.threshold = 3
        self.loc = get_test_files()
        self.expected_informative_sites = numpy.load(os.path.join(self.loc,'chr1_918-test-cutoff-values.npy'))

    def test_informativeness_cutoff(self):
        alignment = os.path.join(self.loc, 'informativeness_cutoff.nex')
        observed_informative_sites = get_informative_sites(alignment, self.threshold)
        small_expected_informative_sites= numpy.array([numpy.nan, numpy.nan, 1., 1.])
        assert observed_informative_sites.all() == small_expected_informative_sites.all()

    def test_full_informativeness_cutoff(self):
        alignment = os.path.join(self.loc, 'chr1_918.nex')
        observed_informative_sites = get_informative_sites(alignment, self.threshold)
        assert observed_informative_sites.all() == \
            self.expected_informative_sites.all()

    def test_cull_informative_rates(self):
        alignment = os.path.join(self.loc, 'chr1_918.nex')
        expected_culled = numpy.load(os.path.join(self.loc,
            'test-culled-rates.npy'))
        # trim to 100 in length to fit w/ fake uniform data
        observed_informative_sites = get_informative_sites(alignment,
                self.threshold)[:100]
        rates = parse_site_rates(os.path.join(self.loc,
            'test-uniform-draw-weights.rates.json'), 10., test = True)
        observed_culled = cull_uninformative_rates(observed_informative_sites, rates)
        assert expected_culled.all() == observed_culled.all() 

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
