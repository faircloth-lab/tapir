"""
File: db.py
Author: Brant Faircloth

Created by Brant Faircloth on 22 October 2011 18:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sqlite3

def create_probe_db(db_name):
    """docstring for create_probe_database"""
    conn = sqlite3.connect(db_name)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    try:
        c.execute("CREATE TABLE loci (id INTEGER PRIMARY KEY AUTOINCREMENT, locus TEXT)")
        c.execute('''CREATE TABLE net (id INT, time INT, pi FLOAT,
            FOREIGN KEY(id) REFERENCES loci(id) DEFERRABLE INITIALLY
            DEFERRED)''')
        c.execute('''CREATE TABLE discrete (id INT, time INT, pi FLOAT, 
            FOREIGN KEY(id) REFERENCES loci(id) DEFERRABLE INITIALLY
            DEFERRED)''')
        c.execute('''CREATE TABLE interval (id INT, interval TEXT, pi FLOAT,
            error FLOAT, FOREIGN KEY(id) REFERENCES loci(id) DEFERRABLE
            INITIALLY DEFERRED)''')
    except sqlite3.OperationalError as e:
        if e[0] == 'table loci already exists':
            answer = raw_input("\n\tPI database already exists.  Overwrite [Y/n]? ")
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
        name, rates, mean_rate, pi, pi_net, times, epochs = locus
        c.execute("INSERT INTO loci(locus) VALUES (?)",
                (os.path.splitext(os.path.basename(name))[0],))
        key = c.lastrowid
        #pdb.set_trace()
        for k,v in enumerate(pi_net):
            c.execute("INSERT INTO net VALUES (?,?,?)", (key, k,
                v))
        if times:
            for k,v in times.iteritems():
                c.execute("INSERT INTO discrete VALUES (?,?,?)", (key, k, v))
        if epochs:
            for k,v in epochs.iteritems():
                c.execute("INSERT INTO interval VALUES (?,?,?,?)", (key, k,
                v['sum(integral)'], v['sum(error)']))
    return
