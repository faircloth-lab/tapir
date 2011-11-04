#!/usr/bin/env python

import os
import re
import collections
import argparse
import sqlite3

import numpy
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages

from tapir.base import FullPaths

import pdb

# lifted from Coding Horror blog
def natural_sort(l):
    """Sorts an iterable using natural (i.e., human) sort"""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('(\d+)', key)]
    return sorted(l, key=alphanum_key)

def get_colormap(string):
    """Ensures that a given color map is valid in matplotlib and returns it."""
    try:
        colormap = getattr(cm, string)
    except AttributeError:
        raise argparse.ArgumentTypeError("Color map not present in matplotlib")
    return colormap

def get_args():
    """Get arguments / options"""
    parser = argparse.ArgumentParser(description="""creates a bar plot of
        phylogenetic informativeness grouped by locus and epoch based on the
        output from pd-ht""")
    parser.add_argument('database', help="SQLite database with pd-ht output",
        action=FullPaths)
    parser.add_argument('--output', help="name of the output file. Format "
        + "will be automatically determined based on the extension.",
        default="pi_by_loci_epoch.pdf")
    parser.add_argument('--font', help="font for plot text. Can be either a "
        + "family (sans-serif) or a specific font name (Helvetica)",
        default="Helvetica")
    parser.add_argument('--colormap', help="matplotlib color map to use",
        default="summer", type=get_colormap)
    parser.add_argument('--width', help="figure width, in inches", default=8,
        type=float)
    parser.add_argument('--height', help="figure height, in inches", default=6,
        type=float)
    parser.add_argument('--keep-extension', help="keep the .nex extension in loci names", action="store_true")
    return parser.parse_args()

def nice_grid(ax):
    """Puts a nice grid behind a plot for the major y axis."""
    ax.yaxis.grid(True, linestyle="-", which="major", color="grey", alpha=0.5)
    ax.set_axisbelow(True)

def create_rects_plot(data, fig, colormap):
    """Creates a bar plot of the PI grouped by locus and epoch"""
    # Each epoch has a width of "1", including an empty spacer bar.
    width = 1 / float(len(data.keys()) + 1)

    # adjust axes to leave room for legend
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.8])

    rects = []  # stores each bar, necessary for correct legend labels
    counter = 0 # tracks colors, necessary for correct tick labels
    for loci, epoch in data.iteritems():
        # Picks a color for this loci
        color = colormap(counter / float(len(data.keys())), 1)

        # Ensure we plot epochs in a correct and consistent order
        epochs_sorted = natural_sort(epoch.keys())
        epoch_values = [epoch[x] for x in epochs_sorted]

        index = numpy.arange(len(epoch))
        rects.append(ax.bar(index + counter * width, epoch_values, width, color=color))
        counter += 1

    # Make a nice grid behind the plots
    nice_grid(ax)

    ax.set_ylabel("Phylogenetic informativeness")
    ax.set_title("Phylogenetic informativeness by loci and epoch")

    ax.set_xticks(index + width * counter / 2)
    ax.set_xticklabels(epochs_sorted)
    ax.set_xlabel("Epochs")

    # Place legend outside the axes
    ax.legend(rects, data.keys(), loc=(1.03, 0.0), title="Loci")

def create_box_plot(data, fig, colors):
    """Creates a box plot by epoch"""
    epochs = []
    loci = []
    pi = []
    tmp = collections.defaultdict(list)
    for locus, epoch in data.iteritems():
        loci.append(locus)
        for key, value in epoch.iteritems():
            tmp[key].append(value)

    for key in natural_sort(tmp.keys()):
        epochs.append(key)
        pi.append(tmp[key])

    ax = fig.add_axes([0.1, 0.1, 0.7, 0.8])
    boxes = ax.boxplot(pi, notch=0, sym='', vert=1, whis=1.5)
    # get the location of the vertical center of each boxplot
    medians = [numpy.average(x.get_xdata()) for x in boxes["medians"]]
    # overplot with points
    # TODO: make colors not give seizure
    points = ax.plot(medians, pi, marker="o", ls='None')

    # make legend
    ax.legend(points, loci, loc=(1.03, 0.0))

    # prettiness etc
    ax.set_xticklabels(epochs)
    ax.set_xlabel("Epochs")
    ax.set_ylabel("Phylogenetic informativeness")
    ax.set_title("PI by loci and epoch")
    nice_grid(ax)

def create_line_plot(data, fig, colors):
    ax = fig.add_axes([0.1, 0.1, 0.7, 0.8])
    all_pi = []
    for locus, pi in data.iteritems():
        x = []
        for a, b in sorted(pi.iteritems()):
            x.append(b)
        all_pi.append(x)

    # transpose because plot() processes 2d array columnwise
    all_pi = numpy.transpose(all_pi)
    points = ax.plot(all_pi)

    ax.set_xlim(ax.get_xlim()[::-1]) # reverse x axis
    ax.legend(points, data.keys(), loc=(1.03, 0.0))
    ax.set_xlabel("MYA")
    ax.set_ylabel("net PI")
    ax.set_title("Net PI over time")
    nice_grid(ax)


def get_epochs(conn):
    c = conn.cursor()
    c.execute("""SELECT loci.locus, interval, pi
                  FROM interval JOIN loci USING (id)""")

    # Rearrange data into the following format:
    # results = {loci1: {epoch1: PI, epoch2: PI, ...}, ...}
    results = collections.defaultdict(dict)
    for row in c:
        loci, epoch, pi = list(row)
        results[loci][epoch] = pi
    return results

def get_net_pi(conn):
    c = conn.cursor()
    c.execute("""SELECT loci.locus, time, pi
                 FROM net JOIN loci USING (id)""")

    results = collections.defaultdict(dict)
    for row in c:
        loci, mya, pi = list(row)
        results[loci][mya] = pi
    return results

def rm_ext(dd):
    """Removes the ".nex" extension from keys in the dictionary `dd`."""
    keys = dd.keys()
    for key in keys:
        new = key.split('.')[0]
        dd[new] = dd.pop(key)
    return dd

def main():
    args = get_args()
    conn = sqlite3.connect(args.database)
    epoch_data = get_epochs(conn)
    net_pi_data = get_net_pi(conn)

    if not args.keep_extension:
        epoch_data, net_pi_data = [rm_ext(x) for x in [epoch_data, net_pi_data]]

    # set defaults
    matplotlib.rc('font', family=args.font)
    matplotlib.rc('figure', figsize=[args.width, args.height])

    # make plots
    create_rects_plot(epoch_data, pyplot.figure(), args.colormap)
    create_box_plot(epoch_data, pyplot.figure(), args.colormap)
    create_line_plot(net_pi_data, pyplot.figure(), args.colormap)

    # For each plot, create a separate image (or separate page for PDFs)
    # TODO: Make this work for more than PDFs
    pdf = PdfPages(args.output)
    for num in pyplot.get_fignums():
        pyplot.figure(num)
        pyplot.savefig(pdf, format="pdf")
    pdf.close()

if __name__ == "__main__":
    main()
