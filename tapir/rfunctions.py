"""
File: rfunctions.py
Author: Brant Faircloth

Created by Brant Faircloth on 22 October 2011 11:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: support functions for R plotting in picme

"""

import rpy2.robjects as robjects
from rpy2.rinterface import RRuntimeError
from rpy2.robjects.packages import importr

import pdb

def load_ggplot(e):
    """if ggplot2 library not present, install"""
    if '''no item called "package:ggplot2"''' in e[0]:
        answer = raw_input("\nggplot2 is not installed.  Install it (and "+
                "dependencies) [Y/n]? ")
        if answer == "Y" or "YES":
            robjects.r('''install.packages("ggplot2")''')
            # there is an rpy2 interface to ggplot2, so let's use it
            import rpy2.robjects.lib.ggplot2 as ggplot2
        else:
            sys.exit(2)
    else:
        print e[0]


def load_sqlite():
    """load the RSQLite library.  install if not present"""
    try:
        robjects.r('library(RSQLite)')
    except RRuntimeError as e:
        if "there is no package called" in e[0]:
            answer = raw_input("\nRSQLite is not installed.  Install it [Y/n]? ")
            if answer == "Y" or "YES":
                # mirror selection choices should pass through on the CLI
                robjects.r('''install.packages("RSQLite")''')
                # load up library after installing:w
                load_sqlite()
            else:
                sys.exit(2)
        else:
            print e[0]

def get_db_conn(db, count = ''):
    """open a connection to an sqlite db in R using SQLite"""
    robjects.r('''drv <- dbDriver("SQLite")''')
    conn_string = '''con{} <- dbConnect(drv, dbname = "{}")'''.format(count, db)
    robjects.r(conn_string)

def close_db_conn():
    """close the db conn"""
    robjects.r('''dbDisconnect(con)''')

def remove_nexus_from_name(frame, col):
    """sometimes locus names contain nexus, remove that"""
    for k,v in enumerate(frame[col]):
        #pdb.set_trace()
        frame[col][k] = v.strip('.nexus')
    return frame

def setup_plotter(name, typ = "png", width = 4, height = 4, units = "in", 
        ppi = 300):
    """we need to setup plotting to the output device (a file, in our case)"""
    grdevices = importr('grDevices')
    if typ == 'png':
        grdevices.png(file = name, width = width, height = height, units=units,
                res = ppi)
        return grdevices
    elif typ == 'pdf':
        grdevices.pdf(file = name, width = width, height = height)
        return grdevices
    elif typ == 'tiff':
        grdevices.tiff(file = name, width = width, height = height, units =
                units, res = ppi, type = "quartz", antialias = "default")
        return grdevices
    elif typ == 'jpeg' or typ == 'jpg':
        grdevices.jpeg(file = name, width = width, height = height, 
                units = units, res=ppi)
        return grdevices
    else:
        raise TypeError("Unknown output format `{0}`".format(typ))
