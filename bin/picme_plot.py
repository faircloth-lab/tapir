"""
File: create_r_plots.py
Author: Brant Faircloth

Created by Brant Faircloth on 22 October 2011 08:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: create plots from data in SQLite

"""

from picme import *
from picme.rfunctions import *

from rpy2 import robjects
from rpy2.rinterface import RRuntimeError

try:
    import rpy2.robjects.lib.ggplot2 as ggplot2
except RRuntimeError as e:
    load_ggplot(e)

import pdb

# set table names as globals, because I'm likely to change them
PI = 'net_informativeness'
LOCUS = 'loci'
DISCRETE = 'time'
INTERVAL = 'epoch'

def get_args():
    """Get arguments / options"""
    parser = argparse.ArgumentParser(description="""Creates plots of
        phylogenetic informativeness based on the
        output of picme.py""")
    parser.add_argument('database', help="SQLite database containing picme output",
        action=FullPaths)
    parser.add_argument('plot_type',help="The type of plot to make")
    # single-locus-net-pi
    # multiple-locus-net-pi
    # 
    parser.add_argument('--loci', help = "", type=get_strings_from_items)
    parser.add_argument('--epoch',help = "", type=get_list_from_ranges)

    parser.add_argument('--output', help="Name of the output file. Format "
        + "will be automatically determined based on the extension. Format "
        + "choices include PDF, PNG, and TIFF",
        default="pi.png")
    parser.add_argument('--width', help="Figure width, in inches", default=8,
        type=float)
    parser.add_argument('--height', help="Figure height, in inches", default=6,
        type=float)
    parser.add_argument('--dpi', help="Figure dpi", default=300,
        type=int)

    return parser.parse_args()

def single_locus_net_informativeness(locus_table, net_pi_table, locus):
    qry = '''"SELECT {0}.locus, mya, pi FROM {0}, {1} 
    WHERE {0}.id = {1}.id AND locus like '{2}'"'''.format(locus_table,
            net_pi_table, locus)
    frame = robjects.r('''dbGetQuery(con, {})'''.format(qry))
    gg_frame = ggplot2.ggplot(frame)
    plot = gg_frame + ggplot2.aes_string(x = 'mya', y='pi') + \
            ggplot2.geom_point(size = 3, alpha = 0.4) + \
            ggplot2.scale_x_reverse() + ggplot2.opts(title = locus)
    return plot

def multiple_locus_net_informativeness_scatterplot(locus_table, net_pi_table,
        loci):
    if loci[0].lower() != 'all':
        qry = '''"SELECT {0}.locus, mya, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id and locus in {2}"'''.format(locus_table,
            net_pi_table, tuple(loci))
    else:
        qry = '''"SELECT {0}.locus, mya, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id"'''.format(locus_table,
            net_pi_table)
    frame = robjects.r('''dbGetQuery(con, {})'''.format(qry))
    gg_frame = ggplot2.ggplot(frame)
    plot = gg_frame + ggplot2.aes_string(x = 'mya', y = 'pi') + \
            ggplot2.geom_point(ggplot2.aes_string(colour = 'locus'), \
            size = 3, alpha = 0.4) + ggplot2.scale_x_reverse()
    return plot


def multiple_locus_net_informativeness_facet(locus_table, net_pi_table, loci):
    if loci[0].lower() != 'all':
        qry = '''"SELECT {0}.locus, mya, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id and locus in {2}"'''.format(locus_table,
            net_pi_table, tuple(loci))
    else:
        qry = '''"SELECT {0}.locus, mya, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id"'''.format(locus_table,
            net_pi_table)
    frame = robjects.r('''dbGetQuery(con, {})'''.format(qry))
    gg_frame = ggplot2.ggplot(frame)
    plot = gg_frame + ggplot2.aes_string(x = 'mya', y='pi') + \
        ggplot2.geom_point(ggplot2.aes_string(colour = 'locus'), size = 3, \
        alpha = 0.4) + ggplot2.scale_x_reverse() + \
        ggplot2.facet_wrap(robjects.Formula('~ locus')) + \
        + ggplot2.opts(**{'legend.position' : 'none'})
    return plot

def make_plot(args):
    plots = []
    if args.plot_type == 'single-locus-net-pi' and args.loci is not None:
        for locus in args.loci:
            plots.append(single_locus_net_informativeness(LOCUS, PI, locus))
    elif args.plot_type == 'multiple-locus-net-pi-facet' and \
            args.loci is not None:
        plots.append(multiple_locus_net_informativeness_facet(LOCUS, PI,
            args.loci))
    elif args.plot_type == 'multiple-locus-net-pi-scatterplot' and \
            args.loci is not None:
        plots.append(multiple_locus_net_informativeness_scatterplot(LOCUS, PI,
            args.loci))

    plotter = setup_plotter(args.output, get_output_type(args.output), args.width,
            args.height, args.dpi)
    for plot in plots:
        plot.plot()
    plotter.dev_off()

def main():
    args = get_args()
    # ggplot2 gets loaded as a module.  here load sqlite
    # and iterface through robjects
    load_sqlite()
    # connect R to db
    get_db_conn(args.database)
    make_plot(args)
    close_db_conn()

if __name__ == '__main__':
    main()
