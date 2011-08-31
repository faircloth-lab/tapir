#!/usr/bin/env python

import argparse
import math
from Bio import Phylo

import pdb

def get_args():
    """ gets command line arguments """
    class T: pass
    a = T()
    a.tree = "test/test-data/Euteleost.tree"
    a.output = "corrected.tree"
    return a

def main():
    args = get_args()
    accumulator = []
    for tree in Phylo.parse(args.tree, 'nexus'):
        tree_length = tree.total_branch_length()
        total_branches = 2 * tree.count_terminals() - 3
        mean_branch_length = tree_length / total_branches
        string_length = math.ceil(math.log10(round(mean_branch_length)))
        if string_length > 1:
            correction_factor = math.pow(10, string_length)
        else:
            correction_factor = 1
        for leaf in tree.get_terminals():
            leaf.branch_length /= correction_factor
        accumulator.append(tree)
    Phylo.write(accumulator, args.output, 'nexus') 


if __name__ == '__main__':
    main()
