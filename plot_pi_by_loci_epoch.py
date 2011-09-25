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

# from estimate_p_i.py, possible TODO: create common library
class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

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
    return parser.parse_args()

def main():
    args = get_args()
    conn = sqlite3.connect(args.database)
    c = conn.cursor()
    c.execute("""SELECT loci.locus, epoch.epoch, epoch.sum_integral
                 FROM loci, epoch
                 WHERE loci.id = epoch.id""")

    # Rearrange data into the following format:
    # results = {loci1: {epoch1: PI, epoch2: PI, ...}, ...}
    results = collections.defaultdict(dict)
    for row in c:
        # 0 = loci; 1 = epoch; 2 = PI
        results[row[0]][row[1]] = row[2]

    # set fonts
    matplotlib.rc('font', family=args.font)

    # Each epoch has a width of "1", including an empty spacer bar.
    width = 1 / float(len(results.keys()) + 1)

    # Plot in the left box on a 2x1 grid to leave space for the legend
    ax = pyplot.figure().add_subplot(121)

    rects = []  # stores each bar, necessary for correct legend labels
    counter = 0 # tracks colors, necessary for correct tick labels
    for loci, epoch in results.iteritems():
        # Picks a color for this loci
        color = args.colormap(counter / float(len(results.keys())), 1)

        # Ensure we plot epochs in a correct and consistent order
        epochs_sorted = natural_sort(epoch.keys())
        epoch_values = [epoch[x] for x in epochs_sorted]

        index = numpy.arange(len(epoch))
        rects.append(ax.bar(index + counter * width, epoch_values, width, color=color))
        counter += 1

    # increase the width of the plot because the legend does not use all
    # of the extra space allocated to it on the 2x1 grid
    pos = list(ax.get_position().bounds)
    pos[2] *= 1.5
    ax.set_position(pos)

    ax.set_ylabel("Phylogenetic informativeness")
    ax.set_title("Phylogenetic informativeness by loci and epoch")

    # XXX: unsure about the math but it seems to work fine
    ax.set_xticks(index + width * counter / 2)
    ax.set_xticklabels(epochs_sorted)
    ax.set_xlabel("Epochs")

    # Place legend at the topright outside corner with 0.05 units of padding
    ax.legend(rects, results.keys(), loc=2, title="Loci",
        bbox_to_anchor=(1.05, 1))

    pyplot.savefig(args.output)

if __name__ == "__main__":
    main()




