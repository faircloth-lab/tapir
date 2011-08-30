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
        self.rates = parse_site_rates('test/test-data/test-uniform-draw-weights.csv')
    
    def runTest(self):
        assert self.rates.all() == numpy.load('test/test-data/test-parsed-rates.npy').all()
    
    def tearDown(self):
        pass

class TestTownsendAndIntegration(unittest.TestCase):
    def setUp(self):
        self.rates = parse_site_rates('test/test-data/test-uniform-draw-weights.csv')
        self.time = numpy.array(range(0, 100 + 1, 1))
    
    def test_townsend_computation(self):
        townsend = get_townsend_pi(self.time, numpy.array([self.rates,]*len(self.time)))
        assert townsend.all() == numpy.load('test/test-data/test-townsend.npy').all()
    
    def test_integration(self):
        vec_integrate = vectorize(integrate)
        # scipy.integrate returns tuple of (integral, upper-error-bound)
        integral, error = vec_integrate(30, 50, self.rates)
        assert integral.all() == numpy.load('test/test-data/test-30-50-integral.npy').all()
        
    def tearDown(self):
        pass
        
if __name__ == '__main__':
    unittest.main()