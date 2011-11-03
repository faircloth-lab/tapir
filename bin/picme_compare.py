#!/usr/bin/env python

"""
File: picme_compare.py
Author: Brant Faircloth

Created by Brant Faircloth on 01 November 2011 14:11 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: Compare two data sets produced

"""

import os
import argparse
import sqlite3

import picme
import picme.rfunctions

from rpy2 import robjects
from rpy2.rinterface import RRuntimeError

try:
    import rpy2.robjects.lib.ggplot2 as ggplot2
except RRuntimeError as e:
    picme.rfunctions.load_ggplot(e)

import pdb

# set table names as globals, because I'm likely to change them
PI = 'net'
LOCUS = 'loci'
DISCRETE = 'discrete'
INTERVAL = 'interval'

def get_args():
    """Get arguments / options"""
    parser = argparse.ArgumentParser(description="""Creates plots comparing
            phylogenetic informativeness of two databases
            output by picme_compute.py""")
    parser.add_argument('db', help="SQLite database containing picme output",
        type=picme.to_full_paths, nargs="+")
    parser.add_argument('plot_type',help="The type of plot to make")
    parser.add_argument('--top', help="""The number of best results across which
            to average""", default = '10')
    parser.add_argument('--db-names', help="SQLite database containing picme output", default=[], type=picme.get_strings_from_items)
    parser.add_argument('--loci', help = "The loci to use in the plot", 
        type=picme.get_strings_from_items, default = None)
    parser.add_argument('--intervals', help = "The intervals to use in the "
        "+ plot.  Must be same between both databases.", 
        type=picme.get_strings_from_items, default = None)
    parser.add_argument('--output', help="Name of the output file. Format "
        + "will be automatically determined based on the extension. Format "
        + "choices include PDF, PNG, and TIFF",
        default=os.path.join(os.getcwd(), "pi.png"))
    parser.add_argument('--width', help="Figure width, in inches",
            default = 8.0, type=float)
    parser.add_argument('--height', help="Figure height, in inches",
            default = 6.0, type=float)
    parser.add_argument('--dpi', help="Figure dpi", default = 150,
        type=int)
    args = parser.parse_args()
    for idx, db in enumerate(args.db):
        try:
            args.db_names[idx]
        except IndexError:
            args.db_names.append(db)
    return args

def order_intervals(intervals):
    e = {int(e.split('-')[0]):e for e in intervals}
    ky = e.keys()
    ky.sort()
    sorted_intervals = [e[k] for k in ky]
    return "c{0}".format(tuple(sorted_intervals))

def get_interval_query(interval, locus_table, interval_table, name, rows):
    if rows == 'all':
        qry = '''"SELECT {0}.locus, interval, pi, '{3}' as db FROM {0}, {1} 
            WHERE {0}.id = {1}.id and interval='{2}' ORDER BY 
            pi DESC"'''.format(locus_table,
            interval_table, interval, name, rows)
    else:
        qry = '''"SELECT {0}.locus, interval, pi, '{3}' as db FROM {0}, {1} 
            WHERE {0}.id = {1}.id and interval='{2}' ORDER BY 
            pi DESC limit {4}"'''.format(locus_table,
            interval_table, interval, name, rows)
    return qry

def get_r_data_by_top(locus_table, interval_table, intervals, names,
            rows):
    # it is faster here to iterate over intervals rather than try and do this
    # with a single query.
    for con,name in names.iteritems():
        for k, interval in enumerate(intervals):
            qry = get_interval_query(interval, locus_table, interval_table, name, rows)
            robjects.r('''con{0}_data{1} <- dbGetQuery(con{0}, {2})'''.format(con, k, qry))
    # prepare string for binding data sets into single stack
    d_string = ','.join(["con{0}_data{1}".format(i,j) for i in range(1, len(names) + 1)
            for j in range(k + 1)])
    #pdb.set_trace()
    frame = robjects.r('''data <- rbind({})'''.format(d_string))
    return frame

def compare_mean_boxplot(locus_table, interval_table, intervals, loci, names, rows):
    frame = get_r_data_by_top(locus_table, interval_table, intervals, names,
            rows)
    if len(intervals) > 1:
        sort_string = '''data$interval <- factor(data$interval, {})'''.format(order_intervals(frame[1]))
        robjects.r(sort_string)
    gg_frame = ggplot2.ggplot(robjects.r('''data'''))
    plot = gg_frame + ggplot2.aes_string(x = 'interval', y = 'pi') + \
                ggplot2.geom_boxplot(ggplot2.aes_string(fill = 'factor(db)'), **{
                    'outlier.size':3,
                    'outlier.colour':'#767676',
                    'outlier.alpha':0.3,
                    'alpha':0.6
                    }
                ) + \
                ggplot2.scale_y_continuous('mean phylogenetic informativeness') + \
                ggplot2.scale_x_discrete('interval (years ago)') + \
                ggplot2.scale_fill_brewer("database", palette='Blues')
    return plot

def compare_sum_barplot(locus_table, interval_table, intervals, loci, names,
        rows):
    frame = get_r_data_by_top(locus_table, interval_table, intervals, names,
            rows)
    #pdb.set_trace()
    frame2 = robjects.r('''agg_data <- aggregate(pi ~ interval + db, data = data, sum)''')
    if len(intervals) > 1:
        sort_string = '''agg_data$interval <- factor(agg_data$interval,{})'''.format(order_intervals(frame2[0]))
        robjects.r(sort_string)
    gg_frame = ggplot2.ggplot(robjects.r('''agg_data'''))
    plot = gg_frame + \
        ggplot2.aes_string(
                x = 'interval', 
                y = 'pi',
                fill='factor(db)'
            ) + \
        ggplot2.geom_bar(**{
            'position':'dodge',
            'colour':'#767676',
            'alpha':0.6
            }
        ) + \
        ggplot2.scale_y_continuous('net phylogenetic informativeness') + \
        ggplot2.scale_x_discrete('interval (years ago)') + \
        ggplot2.scale_fill_brewer("database", palette="Blues")
    return plot

def make_plot(args, names):
    plots = []
    if args.plot_type == 'compare-mean-boxplot':
        plots.append(compare_mean_boxplot(LOCUS, INTERVAL, args.intervals,
            args.loci, names, args.top))
    elif args.plot_type == 'compare-sum-barplot':
        plots.append(compare_sum_barplot(LOCUS, INTERVAL, args.intervals,
            args.loci, names, args.top))
    plotter = picme.rfunctions.setup_plotter(args.output, picme.get_output_type(args.output),
            args.width, args.height, "in", args.dpi)
    for plot in plots:
        plot.plot()
    plotter.dev_off()

def get_or_check_intervals(args, get=False):
    results = []
    for db in args.db:
        con = sqlite3.connect(db)
        cur = con.cursor()
        cur.execute('''SELECT DISTINCT(interval) FROM interval ORDER BY
                interval''')
        results.append(cur.fetchall())
        cur.close()
    assert results[0] == results[1], "Interval values are not equal."
    if get:
        return [str(i[0]) for i in results[0]]

def main():
    # suppress warnings - this is primarily to suppress scale_fill_brewer being
    # unhappy w/ < 3 factors
    robjects.r('''options(warn=-1)''')
    args = get_args()
    if args.intervals == ['all'] or not args.intervals:
        args.intervals = get_or_check_intervals(args, get = True)
    else:
        get_or_check_intervals(args)
    # ggplot2 gets loaded as a module.  here load sqlite
    # and iterface through robjects
    picme.rfunctions.load_sqlite()
    # Adds database connections `con1` through `conn` in the R namespace
    names = {}
    for idx, db in enumerate(args.db, start=1):
        picme.rfunctions.get_db_conn(db, idx)
        names[idx] = args.db_names[idx - 1] # python starts at 0
    make_plot(args, names)

if __name__ == '__main__':
    main()
