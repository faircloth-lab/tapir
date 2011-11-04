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

#import pdb

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def to_full_paths(string):
    return os.path.abspath(os.path.expanduser(string))

def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def mkdir(path):
    try:
        os.mkdir(path)
    except OSError as e:
        pass

def create_unique_dir(path, limit=100):
    """Attempts to create a directory `path`. Returns the name of the
    directory actually created, which may or may not be the same as `path`.

    e.g., if my_directory already exists, it tries to create my_directory.1,
    my_directory.2, ... until my_directory.`limit` is reached.

    Race conditions are possible.
    """
    original = path
    count = 1
    while count < limit:
        try:
            os.mkdir(path)
            return path
        except OSError as e:
            if e.errno == 17: # file exists
                path = "{0}.{1}".format(original, count)
                count += 1
            else:
                raise
    else:
        msg = "could not uniquely create directory {0}: limit `{1}` reached"
        raise Exception(msg.format(original, limit))

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
    except StandardError as e:
        msg = "Cannot convert {0} to a list of integers: {1}"
        raise argparse.ArgumentTypeError(msg.format(name, e))
    return times

def get_strings_from_items(string, name = 'locus'):
    """Convert times input as string to a list"""
    try:
        times = [str(i) for i in string.split(',')]
    except StandardError as e:
        msg = "Cannot convert {0} to a list of loci: {1}"
        raise argparse.ArgumentTypeError(msg.format(name, e))
    return times

def get_list_from_ranges(string):
    """Convert ranges entered as string to nested list"""
    try:
        ranges = [[int(j) for j in i.split('-')] for i in string.split(',')]
    except StandardError as e:
        msg = "Cannot convert spans to a list of integers: {0}"
        raise argparse.ArgumentTypeError(msg.format(e))
    return ranges

def get_files(d, extension):
    if ',' in extension:
        extension = extension.strip(' ').split(',')
    else:
        extension = [extension]
    files = []
    for e in extension:
        files.extend(glob.glob(os.path.join(d, e)))
    if files == []:
        msg = "There appear to be no files of {0} type in {1}"
        raise IOError(msg.format(extension, d))
    else:
        return files
