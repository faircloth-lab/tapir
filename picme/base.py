"""
File: common.py
Author: Brant Faircloth

Created by Brant Faircloth on 22 October 2011 18:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: functions common to picme and helpers

"""

import os
import sys
import glob
import argparse

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def get_output_type(name):
    """get extension from filename"""
    ext = os.path.splitext(name)[1].lstrip('.').lower()
    assert ext in ['pdf','png','tiff','jpeg', 'jpg'], "Filetype must be " + \
        "one of pdf, png, tiff, or jpeg"
    return ext

def get_list_from_ints(string, name = 'time'):
    """Convert times input as string to a list"""
    try:
        times = [int(i) for i in string.split(',')]
    except Exception:
        raise argparse.ArgumentTypeError("Cannot convert {} to list of " + \
                "integers".format(name))
    return times

def get_strings_from_items(string, name = 'locus'):
    """Convert times input as string to a list"""
    try:
        times = [str(i) for i in string.split(',')]
    except Exception:
        raise argparse.ArgumentTypeError("Cannot convert {} to list of " + \
                "loci".format(name))
    return times

def get_list_from_ranges(string):
    """Convert ranges entered as string to nested list"""
    try:
        ranges = [[int(j) for j in i.split('-')] for i in string.split(',')]
    except:
        raise argparse.ArgumentTypeError("Cannot convert spans to list of integers")
    return ranges

def get_files(d, extension):
    if ',' in extension:
        extension = extension.strip(' ').split(',')
    files = []
    for e in extension:
        files.extend(glob.glob(os.path.join(d, e)))
    if files == []:
        print "There appear to be no files of {0} type in {1}".format(extension, d)
        sys.exit(2)
    else:
        return files
