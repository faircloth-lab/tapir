#!/usr/bin/env python

"""
File: tapir_plot.py
Author: Brant Faircloth

Created by Brant Faircloth on 22 October 2011 08:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: create plots from data in SQLite

"""

import os
import argparse

import tapir
import tapir.rfunctions

from rpy2 import robjects
from rpy2.rinterface import RRuntimeError

try:
    import rpy2.robjects.lib.ggplot2 as ggplot2
except RRuntimeError as e:
    tapir.rfunctions.load_ggplot(e)

import pdb

# set table names as globals, because I'm likely to change them
PI = 'net'
LOCUS = 'loci'
DISCRETE = 'discrete'
INTERVAL = 'interval'

def get_args():
    """Get arguments / options"""
    parser = argparse.ArgumentParser(description="""Creates plots of
        phylogenetic informativeness based on the
        output of tapir.py""")
    parser.add_argument('database', help="SQLite database containing tapir output",
        action=tapir.FullPaths)
    parser.add_argument('plot_type',help="The type of plot to make")
    # single-locus-net-pi
    # multiple-locus-net-pi
    # 
    parser.add_argument('--loci', help = "", type=tapir.get_strings_from_items, \
            default = None)
    parser.add_argument('--intervals', help = "", type=tapir.get_strings_from_items, \
            default = None)
    parser.add_argument('--output', help="Name of the output file. Format "
        + "will be automatically determined based on the extension. Format "
        + "choices include PDF, PNG, and TIFF",
        default="pi.png")
    parser.add_argument('--width', help="Figure width, in inches",
            default = 8.0, type=float)
    parser.add_argument('--height', help="Figure height, in inches",
            default = 6.0, type=float)
    parser.add_argument('--dpi', help="Figure dpi", default = 150,
        type=int)
    return parser.parse_args()

def single_locus_net_informativeness(locus_table, net_pi_table, locus):
    qry = '''"SELECT {0}.locus, time, pi FROM {0}, {1} 
    WHERE {0}.id = {1}.id AND locus = '{2}'"'''.format(locus_table,
            net_pi_table, locus)
    frame = robjects.r('''dbGetQuery(con, {})'''.format(qry))
    gg_frame = ggplot2.ggplot(frame)
    plot = gg_frame + ggplot2.aes_string(x = 'time', y='pi') + \
            ggplot2.geom_point(size = 3, alpha = 0.4) + \
            ggplot2.scale_x_reverse('years ago') + \
            ggplot2.scale_y_continuous('phylogenetic informativeness') + \
            ggplot2.opts(title = locus)

    return plot

def multiple_locus_net_informativeness_scatterplot(locus_table, net_pi_table,
        loci):
    if loci[0].lower() != 'all':
        qry = '''"SELECT {0}.locus, time, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id and locus in {2}"'''.format(locus_table,
            net_pi_table, tuple(loci))
    else:
        qry = '''"SELECT {0}.locus, time, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id"'''.format(locus_table,
            net_pi_table)
    frame = robjects.r('''dbGetQuery(con, {})'''.format(qry))
    gg_frame = ggplot2.ggplot(frame)
    plot = gg_frame + ggplot2.aes_string(x = 'time', y = 'pi') + \
            ggplot2.geom_point(ggplot2.aes_string(colour = 'locus'), \
            size = 3, alpha = 0.4) + ggplot2.scale_x_reverse('years ago') + \
            ggplot2.scale_y_continuous('phylogenetic informativeness')
    return plot

def multiple_locus_net_informativeness_facet(locus_table, net_pi_table, loci):
    if loci[0].lower() != 'all':
        qry = '''"SELECT {0}.locus, time, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id and locus in {2}"'''.format(locus_table,
            net_pi_table, tuple(loci))
    else:
        qry = '''"SELECT {0}.locus, time, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id"'''.format(locus_table,
            net_pi_table)
    frame = robjects.r('''dbGetQuery(con, {})'''.format(qry))
    gg_frame = ggplot2.ggplot(frame)
    plot = gg_frame + ggplot2.aes_string(x = 'time', y='pi') + \
        ggplot2.geom_point(ggplot2.aes_string(colour = 'locus'), size = 3, \
        alpha = 0.4) + ggplot2.scale_x_reverse('years ago') + \
        ggplot2.facet_wrap(robjects.Formula('~ locus')) + \
        ggplot2.opts(**{'legend.position' : 'none'}) + \
        ggplot2.scale_y_continuous('phylogenetic informativeness')
    return plot

def order_intervals(intervals):
    e = {int(e.split('-')[0]):e for e in intervals}
    ky = e.keys()
    ky.sort()
    sorted_intervals = [e[k] for k in ky]
    return "c{0}".format(tuple(sorted_intervals))

def get_interval_query(intervals, loci, locus_table, interval_table):
    #pdb.set_trace()
    if intervals[0].lower() != 'all' and (loci is None or loci[0] == 'all'):
        qry = '''"SELECT {0}.locus, interval, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id and interval in {2}"'''.format(locus_table,
            interval_table, tuple(intervals))
    elif intervals[0].lower() != 'all' and loci[0].lower() != 'all':
        qry = '''"SELECT {0}.locus, interval, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id and interval in {2} and locus in {3}"'''.format(locus_table,
            interval_table, tuple(intervals), tuple(loci))
    elif intervals[0].lower() == 'all' and loci[0].lower() != 'all':
        qry = '''"SELECT {0}.locus, interval, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id and locus in {2}"'''.format(locus_table,
            interval_table, tuple(loci))
    else:
        qry = '''"SELECT {0}.locus, interval, pi FROM {0}, {1} 
            WHERE {0}.id = {1}.id"'''.format(locus_table,
            interval_table)
    return qry

def interval(locus_table, interval_table, intervals, loci, boxplot = True):
    qry = get_interval_query(intervals, loci, locus_table, interval_table)
    frame = robjects.r('''data <- dbGetQuery(con, {})'''.format(qry))
    # because we're sorting by interval, which is a factor, we need to
    # explicitly re-sort the data by the first integer value
    # of the interval.  This is a bit cumbersome, because sorting
    # in R is less than pleasant.
    sort_string = '''data$interval <- factor(data$interval, {})'''.format(order_intervals(frame[1]))
    robjects.r(sort_string)
    gg_frame = ggplot2.ggplot(robjects.r('''data'''))
    if boxplot:
        plot = gg_frame + ggplot2.aes_string(x = 'interval', y = 'pi') + \
                ggplot2.geom_boxplot(**{
                    'outlier.size':0, 
                    'alpha':0.3
                    }
                ) + \
                ggplot2.geom_jitter(ggplot2.aes_string(color = 'locus'), size = 3, \
                alpha = 0.6, position=ggplot2.position_jitter(width=0.25)) + \
                ggplot2.scale_y_continuous('phylogenetic informativeness') + \
                ggplot2.scale_x_discrete('interval (years ago)')

    else:
        plot = gg_frame + ggplot2.aes_string(x = 'interval', y = 'pi',
                fill='locus') + ggplot2.geom_bar() + \
                ggplot2.facet_wrap(robjects.Formula('~ locus')) + \
                ggplot2.opts(**{
                    'axis.text.x':ggplot2.theme_text(angle = -90,  hjust = 0),
                    'legend.position':'none'
                    }) + \
                ggplot2.scale_y_continuous('phylogenetic informativeness') + \
                ggplot2.scale_x_discrete('interval (years ago)')
    return plot

def make_plot(args):
    plots = []
    if args.plot_type == 'pi-locus' and args.loci is not None:
        for locus in args.loci:
            plots.append(single_locus_net_informativeness(LOCUS, PI, locus))
    elif args.plot_type == 'pi-facet' and \
            args.loci is not None:
        plots.append(multiple_locus_net_informativeness_facet(LOCUS, PI,
            args.loci))
    elif args.plot_type == 'pi-scatterplot' and \
            args.loci is not None:
        plots.append(multiple_locus_net_informativeness_scatterplot(LOCUS, PI,
            args.loci))
    elif args.plot_type == 'pi-interval-boxplot' and args.intervals is not None:
        plots.append(interval(LOCUS, INTERVAL, args.intervals, args.loci))
    elif args.plot_type == 'pi-interval-barplot' and args.intervals is not None:
        plots.append(interval(LOCUS, INTERVAL, args.intervals, args.loci,
            boxplot = False))

    plotter = tapir.rfunctions.setup_plotter(args.output, tapir.get_output_type(args.output),
            args.width, args.height, "in", args.dpi)
    for plot in plots:
        plot.plot()
    plotter.dev_off()

def main():
    args = get_args()
    # ggplot2 gets loaded as a module.  here load sqlite
    # and iterface through robjects
    tapir.rfunctions.load_sqlite()
    # connect R to db
    tapir.rfunctions.get_db_conn(args.database)
    make_plot(args)
    tapir.rfunctions.close_db_conn()

if __name__ == '__main__':
    main()
