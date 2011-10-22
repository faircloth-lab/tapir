"""
File: rfunctions.py
Author: Brant Faircloth

Created by Brant Faircloth on 22 October 2011 11:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: support functions for R plotting in picme

"""

import rpy2.robjects as robjects
from rpy2.rinterface import RRuntimeError

import pdb

def load_ggplot(e):
    if '''no item called "package:ggplot2"''' in e[0]:
        answer = raw_input("\nggplot2 is not installed.  Install it (and "+
                "dependencies) [Y/n]? ")
        if answer == "Y" or "YES":
            robjects.r('''install.packages("ggplot2", dependencies = TRUE)''')
            # there is an rpy2 interface to ggplot2, so let's use it
            import rpy2.robjects.lib.ggplot2 as ggplot2
        else:
            sys.exit(2)
    else:
        print e[0]


def load_sqlite():
    try:
        robjects.r('library(RSQLite)')
    except RRuntimeError as e:
        if "no package called 'SQLite'" in e[0]:
            answer = raw_input("\nRSQLite is not installed.  Install it [Y/n]? ")
            if answer == "Y" or "YES":
                # mirror selection choices should pass through on the CLI
                robjects.r('''install.packages("RSQLite", dependencies = TRUE)''')
                # load up library after installing:w
                load_sqlite()
            else:
                sys.exit(2)
        else:
            print e[0]

def get_db_conn(db):
    robjects.r('''drv <- dbDriver("SQLite")''')
    conn_string = '''con <- dbConnect(drv, dbname = "{}")'''.format(db)
    robjects.r(conn_string)

def close_db_conn(db):
    pass

def remove_nexus_from_name(frame, col):
    for k,v in enumerate(frame[col]):
        #pdb.set_trace()
        frame[col][k] = v.strip('.nexus')
    return frame
