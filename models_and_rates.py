#!/usr/bin/env python

import os
import sys
import imp
import json
import numpy
import argparse
import dendropy
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
    parser.add_argument('--tree-format', dest='tree_format', help="The format of the tree",
        choices=['nexus','newick'], default='newick')
    parser.add_argument('--output', dest='output', help="The path to the output"
        +" directory", default=os.getcwd(), action=FullPaths)
    parser.add_argument('--hyphy', dest='hyphy', default="hyphy1", help="The "
        +"path to hyphy (if not in $PATH)")
    #parser.add_argument('--test', action='store_true')
    return parser.parse_args()

def parse_site_rates(rate_file, correction = 1):
    """parse the site rate file returned from hyphy to a vector of rates"""
    data = json.load(open(rate_file, 'r'))
    rates = numpy.array([line["rate"] for line in data["sites"]["rates"]])
    return numpy.reshape(rates/correction, (-1,1))

def get_townsend_pi(time, rates):
    return 16 * (rates**2) * time * numpy.exp(-(4 * rates * time))

def integrate(start, stop, rate):
    return quad(get_townsend_pi, start, stop, args=(rate))

def correct_branch_lengths(tree_file, format, d = ""):
    """The Townsend phydesign code corrects branch lengths based on an algorithm
    implement that in python"""
    tree = dendropy.Tree.get_from_path(tree_file, format)
    depth = tree.seed_node.distance_from_tip()
    mean_branch_length = tree.length()/(2 * len(tree.leaf_nodes()) - 3)
    string_len = len(str(int(mean_branch_length + 0.5)))
    if string_len > 1:
        correction_factor = 10 ** string_len
    else:
        correction_factor = 1
    for edge in tree.preorder_edge_iter():
        if edge.length:
            edge.length /= correction_factor
    pth = os.path.join(d, 'Tree_{}_{}.newick'.format(correction_factor, depth))
    tree.write_to_path(pth, 'newick')
    return depth, correction_factor, pth


def main():
    """main loop"""
    args = get_args()
    # correct branch lengths
    depth, correction, tree = correct_branch_lengths(args.tree, args.tree_format, d=args.output)
    # generate a vector of times given start and stops
    time = numpy.array(range(args.start, args.end + 1, args.step))
    hyphy = Popen([args.hyphy, 'templates/models_and_rates.bf'], \
        stdin=PIPE, stdout=PIPE)
    output = os.path.join(args.output, os.path.basename(args.alignment) + '.rates')
    towrite = "\n".join([args.alignment, tree, output])
    #pdb.set_trace()
    stdout, stderr = hyphy.communicate(towrite)
    rates = parse_site_rates(output, correction = correction)
    # send column of times and vector of site rates to get_townsend_pi.
    # Because of structure, we can take advantage of numpy's
    # elementwise speedup
    phylogenetic_informativeness = \
        get_townsend_pi(time, rates)
    # vectorize the integral function to take our rates array as input
    vec_integrate = vectorize(integrate)
    # scipy.integrate returns tuple of (integral, upper-error-bound)
    # TODO:  figure out how we want to handle diff. time intervals here
    integral, error = vec_integrate(0.5, 0.6, rates)
    pdb.set_trace()



if __name__ == '__main__':
    main()