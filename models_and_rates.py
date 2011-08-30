#!/usr/bin/env python

import os
import sys
import argparse
from subprocess import Popen, PIPE

import pdb

def get_args():
    """get CLI arguments and options"""
    parser = argparse.ArgumentParser(description='Select model and generate site rates')
    parser.add_argument('alignment', help="The input alignment")
    parser.add_argument('tree', help="The input tree")
    parser.add_argument('--output', dest='output', help="The path to the output directory", default=os.getcwd())
    parser.add_argument('--hyphy', dest='hyphy', default="HYPHY", help="The path to hyphy (if not in $PATH)")
    #parser.add_argument('--csv', action='store_true')
    return parser.parse_args()

def main():
    """main loop"""
    args = get_args()
    hyphy = Popen([args.hyphy, 'templates/models_and_rates.bf'], stdin=PIPE, stdout=PIPE)
    output = os.path.join(os.path.abspath(os.path.expanduser(args.output)), os.path.basename(args.alignment) + '.rates')
    towrite = "\n".join([os.path.abspath(os.path.expanduser(args.alignment)), os.path.abspath(os.path.expanduser(args.tree)), output])
    #pdb.set_trace()
    stdout, stderr = hyphy.communicate(towrite)
    print 'Done'

if __name__ == '__main__':
    main()