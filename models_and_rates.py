#!/usr/bin/env python

import os
import sys
import time
import glob
import json
import numpy
import sqlite3
import argparse
import dendropy
from collections import defaultdict
from scipy import integrate
from scipy import vectorize
from subprocess import Popen, PIPE


import pdb

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def is_dir(dirname):
    if not os.path.isdir:
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

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
    parser.add_argument('alignments', help="The input alignment",
            action=FullPaths, type = is_dir)
    parser.add_argument('tree', help="The input tree", action=FullPaths)
    parser.add_argument('--times', help="The start time of interest (MYA)", type=get_times_from_list)
    parser.add_argument('--epochs', help="The start time of interest (MYA)", type=get_epochs_from_list)
    parser.add_argument('--tree-format', dest='tree_format', help="The format of the tree",
        choices=['nexus','newick'], default='newick')
    parser.add_argument('--output', dest='output', help="The path to the output"
        +" directory", default=os.getcwd(), action=FullPaths)
    parser.add_argument('--hyphy', dest='hyphy', default="hyphy1", help="The "
        +"path to hyphy (if not in $PATH)")
    parser.add_argument('--threshold', default=3, help="Minimum number of taxa"
        +" without a gap for a site to be considered informative")
    parser.add_argument('--multiprocessing', default = False, action =
        'store_true')
    parser.add_argument('--site-rates', default = False, action = 'store_true')
    #parser.add_argument('--test', action='store_true')
    return parser.parse_args()

def parse_site_rates(rate_file, correction = 1, test = False):
    """Parse the site rate file returned from hyphy to a vector of rates"""
    data = json.load(open(rate_file, 'r'))
    rates = numpy.array([line["rate"] for line in data["sites"]["rates"]])
    corrected = rates/correction
    if not test:
        data["sites"]["corrected_rates"] = [{"site":k + 1,"rate":v} \
                for k,v in enumerate(corrected)]
        json.dump(data, open(rate_file,'w'), indent = 4)
    return corrected

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

def get_informative_sites(alignment, threshold=4):
    """Returns a list, where True indicates a site which was over the threshold
    for informativeness.
    """
    taxa = dendropy.DnaCharacterMatrix.get_from_path(alignment, 'nexus')
    results = defaultdict(int)
    for cells in taxa.vectors():
        assert len(cells) == taxa.vector_size # should all have equal lengths
        for idx, cell in enumerate(cells):
            results[idx] += 1 if str(cell).upper() in "ATGC" else 0
    return numpy.array([results[x] >= threshold for x in sorted(results)])

def cull_uninformative_rates(rates, inform):
    """Zeroes out rates which are uninformative"""
    return rates * inform

def worker(params):
    #pdb.set_trace()
    time_vector, hyphy, template, towrite, output, correction, alignment, times, epochs, threshold = params
    # if twowrite is set, run hyphy, else, we've sent site rates
    sys.stdout.write(".")
    sys.stdout.flush()
    if towrite:
        hyphy = Popen([hyphy, template], stdin=PIPE, stdout=PIPE)
        stdout, stderr = hyphy.communicate(towrite)
        rates = parse_site_rates(output, correction = correction)
        good_sites = get_informative_sites(alignment, threshold)
        rates = cull_uninformative_rates(rates, good_sites)
    else:
        rates = parse_site_rates(output, correction = correction)
    # send column of times and vector of site rates to get_townsend_pi.
    # Because of structure, we can take advantage of numpy's
    # elementwise speedup
    phylogenetic_informativeness = get_townsend_pi(time_vector, rates)
    pi_times = get_net_pi_for_periods(phylogenetic_informativeness, times)
    pi_epochs = get_net_integral_for_epochs(rates, epochs)
    return alignment, pi_times, pi_epochs

def create_probe_db(db_name):
    """docstring for create_probe_database"""
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    try:
        c.execute("CREATE TABLE loci (id INTEGER PRIMARY KEY AUTOINCREMENT, locus TEXT)")
        c.execute('''CREATE TABLE time (id INT, time INT, pi FLOAT, 
            FOREIGN KEY(id) REFERENCES loci(id) DEFERRABLE INITIALLY
            DEFERRED)''')
        c.execute('''CREATE TABLE epoch (id INT, epoch TEXT, sum_integral FLOAT,
            sum_error FLOAT, FOREIGN KEY(id) REFERENCES loci(id) DEFERRABLE
            INITIALLY DEFERRED)''')
    except sqlite3.OperationalError, e:
        if e[0] == 'table loci already exists':
            answer = raw_input("\nPI database already exists.  Overwrite [Y/n]? ")
            if answer == "Y" or "YES":
                os.remove(db_name)
                conn, c = create_probe_db(db_name)
            else:
                sys.exit(2)
        else:
            raise sqlite3.OperationalError 
            pdb.set_trace()
    return conn, c

def insert_pi_data(conn, c, pis):
    for locus in pis:
        name, times, epochs = locus
        c.execute("INSERT INTO loci(locus) values (?)",
                (os.path.basename(name),))
        key = c.lastrowid
        if times:
            for k,v in times.iteritems():
                c.execute("INSERT INTO time values (?,?,?)", (key, k, v))
        if epochs:
            for k,v in epochs.iteritems():
                c.execute("INSERT INTO epoch VALUES (?,?,?,?)", (key, k,
                v['sum(integral)'], v['sum(error)']))
    return

def get_files(d, extension):
    files = glob.glob(os.path.join(d, extension))
    if files == []:
        print "There appear to be no files of {} type in {}".format(extension, d)
        sys.exit(2)
    else:
        return files

def main():
    """Main loop"""
    args = get_args()
    # correct branch lengths
    tree_depth, correction, tree = correct_branch_lengths(args.tree, args.tree_format, d = args.output)
    # generate a vector of times given start and stops
    time_vector = get_time(0, int(tree_depth))
    params = []
    # get path to batch/template file for hyphy
    template = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
            'templates/models_and_rates.bf')
    if not args.site_rates:
        print "Estimating site rates and PI for files:"
        for alignment in get_files(args.alignments, '*.nex'):
            output = os.path.join(args.output, os.path.basename(alignment) + '.rates')
            towrite = "\n".join([alignment, tree, output])
            params.append([time_vector, args.hyphy, template, towrite, output, correction, alignment,
                args.times, args.epochs, args.threshold])
    else:
        print "Estimating PI for files (--site-rate option):"
        for rate_file in get_files(args.alignments, '*.rates'):
            params.append([time_vector, args.hyphy, template, None, rate_file,
                correction, rate_file, args.times, args.epochs,
                args.threshold])
    if not args.multiprocessing:
        pis = map(worker, params)
    else:
        from multiprocessing import Pool, cpu_count
        pool = Pool(processes = cpu_count() - 1)
        pis = pool.map(worker, params)
    # store results somewhere
    db_name = os.path.join(args.output,
        'phyogenetic-informativeness.sqlite')
    conn, c = create_probe_db(db_name)
    insert_pi_data(conn, c, pis)
    conn.commit()
    c.close()
    conn.close()

if __name__ == '__main__':
    main()
