#!/usr/bin/env python

import os
import sys
import imp
import json
import numpy
import argparse
import dendropy
from scipy import integrate
from scipy import vectorize
from subprocess import Popen, PIPE


import pdb

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def get_epochs_from_list(string):
    """Convert epochs/spans entered as string to nested list"""
    try:
        epochs = [[int(j) for j in i.split('-')] for i in string.split(',')]
    except:
        raise argparse.ArgumentTypeError("Cannot convert spans to list of integers")
    return epochs

def get_times_from_list(string):
    """Convert times input as string to a list"""
    try:
        times = [int(i) for i in string.split(',')]
    except:
        raise argparse.ArgumentTypeError("Cannot convert time to list of integers")
    return times

def get_args():
    """Get CLI arguments and options"""
    parser = argparse.ArgumentParser(description='Select model and generate site rates')
    parser.add_argument('alignment', help="The input alignment", action=FullPaths)
    parser.add_argument('tree', help="The input tree", action=FullPaths)
    parser.add_argument('--times', help="The start time of interest (MYA)", type=get_times_from_list)
    parser.add_argument('--epochs', help="The start time of interest (MYA)", type=get_epochs_from_list)
    parser.add_argument('--tree-format', dest='tree_format', help="The format of the tree",
        choices=['nexus','newick'], default='newick')
    parser.add_argument('--output', dest='output', help="The path to the output"
        +" directory", default=os.getcwd(), action=FullPaths)
    parser.add_argument('--hyphy', dest='hyphy', default="hyphy1", help="The "
        +"path to hyphy (if not in $PATH)")
    #parser.add_argument('--test', action='store_true')
    return parser.parse_args()

def parse_site_rates(rate_file, correction = 1):
    """Parse the site rate file returned from hyphy to a vector of rates"""
    data = json.load(open(rate_file, 'r'))
    rates = numpy.array([line["rate"] for line in data["sites"]["rates"]])
    return rates/correction

def get_townsend_pi(time, rates):
    """Townsend et al. Equation10 """
    return 16 * (rates**2) * time * numpy.exp(-(4 * rates * time))

def get_integral_over_times(start, stop, rate):
    """Integrate over start and stop times"""
    return integrate.quad(get_townsend_pi, start, stop, args=(rate))

def get_time(start, stop, step = 1):
    """Given start and stop times, return a column of times over which we're working"""
    # reshape array into columns from row
    return numpy.reshape(numpy.array(range(start, stop, step)), (-1,1))

def correct_branch_lengths(tree_file, format, d = ""):
    """Scale branch lengths to values shorter than 100"""
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

def get_net_pi_for_periods(pi, times):
    """Sum across the PI values for the requested times"""
    sums = numpy.sum(pi, axis=1)[times]
    return dict(zip(times, sums))

def get_net_integral_for_epochs(rates, epochs):
    """Given a set of epochs, integrate rates over those start and stop times"""
    # vectorize the integral function to take our rates array as input
    vec_integrate = vectorize(get_integral_over_times)
    # scipy.integrate returns tuple of (integral, upper-error-bound)
    # TODO:  figure out how we want to handle diff. time intervals here
    epochs_results = {}
    for span in epochs:
        name = "{}-{}".format(span[0],span[1])
        assert span[0] < span[1], \
            "Start time [{}] is sooner than end time [{}]".format(span[0],span[1])
        integral, error = vec_integrate(span[0],span[1], rates)
        epochs_results[name] = {'sum(integral)':sum(integral), 'sum(error)':sum(error)}
    return epochs_results

def main():
    """Main loop"""
    args = get_args()
    # correct branch lengths
    tree_depth, correction, tree = correct_branch_lengths(args.tree, args.tree_format, d = args.output)
    # generate a vector of times given start and stops
    time = get_time(0, int(tree_depth))
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
    phylogenetic_informativeness = get_townsend_pi(time, rates)
    if args.times:
        print get_net_pi_for_periods(phylogenetic_informativeness, args.times)
    if args.epochs:
        print get_net_integral_for_epochs(rates, args.epochs)


if __name__ == '__main__':
    main()