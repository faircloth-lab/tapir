#!/usr/bin/env python

import os
import sys
import numpy
import argparse
from scipy import vectorize
from subprocess import Popen, PIPE
from scipy.integrate import quad

import pdb

class FullPaths(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def get_args():
    """get CLI arguments and options"""
    parser = argparse.ArgumentParser(description='Select model and generate site rates')
    parser.add_argument('alignment', help="The input alignment", action=FullPaths)
    parser.add_argument('tree', help="The input tree", action=FullPaths)
    parser.add_argument('start', help="The start time of interest (MYA)", type=int)
    parser.add_argument('end', help="The end time of interest (MYA)", type=int)
    parser.add_argument('--step', dest='step', help="The step distances between"
        +"`start` and `end`", default=1, type=int)
    parser.add_argument('--output', dest='output', help="The path to the output"
        +" directory", default=os.getcwd(), action=FullPaths)
    parser.add_argument('--hyphy', dest='hyphy', default="HYPHY", help="The "
        +"path to hyphy (if not in $PATH)")
    parser.add_argument('--test', action='store_true')
    return parser.parse_args()

def parse_site_rates(rate_file):
    """parse the site rate file returned from hyphy to a vector of rates"""
    rates = numpy.array([float(line.split(',')[2]) 
        for line in open(rate_file, 'rU') if not line.startswith('site')])
    return numpy.reshape(rates, (-1,1))

def get_townsend_pi(time, rates):
    return 16 * (rates**2) * time * numpy.exp(-(4 * rates * time))

def integrate(start, stop, rate):
    return quad(get_townsend_pi, start, stop, args=(rate))

def main():
    """main loop"""
    args = get_args()
    # generate a vector of times given start and stops
    time = numpy.array(range(args.start, args.end + 1, args.step))
    hyphy = Popen([args.hyphy, 'templates/models_and_rates.bf'], stdin=PIPE, stdout=PIPE)
    output = os.path.join(args.output, os.path.basename(args.alignment) + '.rates')
    towrite = "\n".join([args.alignment, args.tree, output])
    stdout, stderr = hyphy.communicate(towrite)
    rates = parse_site_rates(output)
    # send column of times and vector of site rates to get_townsend_pi.
    # Because of structure, we can take advantage of numpy's
    # elementwise speedup
    townsend = get_townsend_pi(time, numpy.array([rates,]*len(time)))
    # vectorize the integral function to take our rates array as input
    vec_integrate = vectorize(integrate)
    # scipy.integrate returns tuple of (integral, upper-error-bound)
    # TODO:  figure out how we want to handle diff. time intervals here
    integral, error = vec_integrate(30, 50, rates)
    pdb.set_trace()
    

if __name__ == '__main__':
    main()