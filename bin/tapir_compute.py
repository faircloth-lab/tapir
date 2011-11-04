#!/usr/bin/env python


import os
import sys
import glob
import numpy
import sqlite3
import argparse

from subprocess import Popen, PIPE

import tapir


import pdb

def get_args():
    """Get CLI arguments and options"""
    parser = argparse.ArgumentParser(description="""tapir:  compute the
            phylogenetic informativeness of DNA loci""")

    parser.add_argument('alignments', help="The folder of alignments",
        action=tapir.FullPaths, type=tapir.is_dir)
    parser.add_argument('tree', help="The input tree", action=tapir.FullPaths)

    required = parser.add_argument_group("required arguments")
    required.add_argument('--times', help="""Comma-separated list of start
        times of interest (MYA)""", type=tapir.get_list_from_ints, required=True)
    required.add_argument('--intervals', help="""Comma-separated list of
        interval ranges of interest (i,e. in MYA)""", type=tapir.get_list_from_ranges, 
        required=True)

    parser.add_argument('--tree-format', help="The format of the tree",
        dest='tree_format', choices=['nexus','newick'], default='newick')
    parser.add_argument('--output', help="The path to the output directory",
        default=os.getcwd(), action=tapir.FullPaths, type=tapir.is_dir)
    parser.add_argument('--hyphy', help="The path to hyphy (if not in $PATH)",
        default="hyphy2")
    parser.add_argument('--template', help="""Path to the hypy
        temlate file if in non-standard location""", default = None,)
    parser.add_argument('--threshold', help="""Minimum number of taxa without
        a gap for a site to be considered informative""", default=3, type=int)
    parser.add_argument('--multiprocessing', help="""Enable parallel
        calculation of rates""", default=False, action='store_true')
    parser.add_argument('--site-rates', default=False, action='store_true',
        help="Use previously calculated site rates")
    return parser.parse_args()

def welcome_message():
    return '''
    ***************************************************
    *                                                 *
    * tapir: high-throughput estimates of             *
    * phylogenetic informativeness                    *
    *                                                 *
    * (c) 2011 Brant Faircloth, Jonathan Chang,       *
    * Mike Alfaro                                     *
    *                                                 *
    * PhyDesign was created by the Townsend Lab       *
    * (http://phydesign.townsend.yale.edu)            *
    *                                                 *
    * When using tapir, please cite:                  *
    *                                                 *
    *   - J.P. Townsend, 2007. Profiling              *
    *     phylogenetic informativeness. Systematic    *
    *     Biology, 56(2), 222-231.                    *
    *                                                 *
    *   - Pond, S.L.K., Frost, S.D.W., and S.V. Muse, *
    *     2005. Hyphy: hypothesis testing using       *
    *     phylogenies. Bioinformatics, 21(5), 676-9.  *
    *                                                 *
    * Many thanks to Francesc Lopez-Giraldez and      *
    * Jeffrey Townsend for providing us with a copy   *
    * of their web-application source code.           *
    *                                                 *
    ***************************************************\n\n'''

def worker(params):
    #pdb.set_trace()
    time_vector, hyphy, template, towrite, output, correction, alignment, times, epochs, threshold = params
    # if twowrite is set, run hyphy, else, we've sent site rates
    sys.stdout.write(".")
    sys.stdout.flush()
    if towrite:
        try:
            hyphy = Popen([hyphy, template], stdin=PIPE, stdout=PIPE)
        except OSError as e:
            if e.errno == 2:
                raise Exception("hyphy not found")
            raise Exception("couldn't communicate with hyphy: {0}".format(e))
        stdout, stderr = hyphy.communicate(towrite)
        if stdout.startswith("Error"):
            raise Exception("hyphy error: {0}".format(stdout))
        rates = tapir.parse_site_rates(output, correction = correction)
        good_sites = tapir.get_informative_sites(alignment, threshold)
        rates = tapir.cull_uninformative_rates(rates, good_sites)
    else:
        rates = tapir.parse_site_rates(output, correction = correction)
    # compute the mean, ensuring we mask the nans.
    mean_rate = numpy.mean(numpy.ma.masked_array(rates, numpy.isnan(rates)))
    # send column of times and vector of site rates to get_townsend_pi.
    # Because of structure, we can take advantage of numpy's
    # elementwise speedup
    pi = tapir.get_townsend_pi(time_vector, rates)
    # len(rates) gives the number of "#Sites" per PhyDesign website
    # numpy.sum(numpy.isnan(rates)) counts undefined sites
    # numpy.sum(numpy.isfinite(rates)) counts #Rates per PhyDesign website
    # pdb.set_trace()
    pi_net = numpy.nansum(pi, axis = 1)
    pi_times = tapir.get_net_pi_for_periods(pi, times)
    # remove the numpy.nan records before computing the integral
    pi_epochs = tapir.get_net_integral_for_epochs(rates[numpy.isfinite(rates)], epochs)
    return alignment, rates, mean_rate, pi, pi_net, pi_times, pi_epochs

def main():
    """Main loop"""
    args = get_args()
    # print message
    print welcome_message()
    # make output dir
    args.output = tapir.create_unique_dir(args.output)
    # correct branch lengths
    tree_depth, correction, tree = tapir.correct_branch_lengths(args.tree, args.tree_format, d = args.output)
    # generate a vector of times given start and stops
    time_vector = tapir.get_time(0, int(tree_depth))
    params = []
    # get path to batch/template file for hyphy
    if not args.template:
        template = tapir.get_hyphy_conf()
    else:
        template = args.template
    if not args.site_rates:
        print "\nEstimating site rates and PI for files:"
        for alignment in tapir.get_files(args.alignments, '*.nex,*.nexus'):
            output = os.path.join(args.output, os.path.basename(alignment) + '.rates')
            towrite = "\n".join([alignment, tree, output])
            params.append([time_vector, args.hyphy, template, towrite, output, correction, alignment,
                args.times, args.intervals, args.threshold])
    else:
        print "Estimating PI for files (--site-rate option):"
        for rate_file in tapir.get_files(args.alignments, '*.rates'):
            params.append([time_vector, args.hyphy, template, None, rate_file,
                correction, rate_file, args.times, args.intervals,
                args.threshold])
    if not args.multiprocessing:
        pis = map(worker, params)
    else:
        from multiprocessing import Pool, cpu_count
        pool = Pool(processes = cpu_count() - 1)
        pis = pool.map(worker, params)
    # store results somewhere
    db_name = os.path.join(args.output,
        'phylogenetic-informativeness.sqlite')
    sys.stdout.write("\nStoring results in {0}...".format(db_name))
    sys.stdout.flush()
    conn, c = tapir.create_probe_db(db_name)
    tapir.insert_pi_data(conn, c, pis)
    conn.commit()
    sys.stdout.write("DONE")
    sys.stdout.flush()
    print "\n"
    c.close()
    conn.close()

if __name__ == '__main__':
    main()
