"""
File: create_r_plots.py
Author: Brant Faircloth

Created by Brant Faircloth on 22 October 2011 08:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: create plots from data in SQLite

"""

from picme.rfunctions import *

import rpy2.robjects as robjects
from rpy2.rinterface import RRuntimeError
from rpy2.robjects.packages import importr

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
    pass

def single_locus_net_informativeness(locus_table, net_pi_table, locus):
    qry = '''"SELECT {0}.locus, mya, pi FROM {0}, {1} 
    WHERE {0}.id = {1}.id AND locus = '{2}'"'''.format(locus_table,
            net_pi_table, locus)
    frame = robjects.r('''dbGetQuery(con, {})'''.format(qry))
    frame = remove_nexus_from_name(frame, 0)
    gg_frame = ggplot2.ggplot(frame)
    plot = gg_frame + ggplot2.aes_string(x = 'mya', y='pi') + \
            ggplot2.geom_point(size = 3, alpha = 0.4) + \
            ggplot2.scale_x_reverse()
    return plot

def setup_plotter(name = 'test.png', typ = 'png', width = 512, height = 512):
    grdevices = importr('grDevices')
    if typ == 'png':
        grdevices.png(file = name, width = width, height = height)
        return grdevices

def main():
    # ggplot2 gets loaded as a module.  here load sqlite
    # and iterface through robjects
    load_sqlite()
    # connect R to db
    get_db_conn('/Users/bcf/Dropbox/Research/alfaro/faircloth-pi/plot-test/phylogenetic-informativeness.sqlite')
    plot = single_locus_net_informativeness(LOCUS, PI, 'adora3.nexus')
    plotter = setup_plotter()
    plot.plot()
    pdb.set_trace()
    plotter.dev_off()

if __name__ == '__main__':
    main()
