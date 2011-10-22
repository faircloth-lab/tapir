"""
File: __init__.py
Author: Brant Faircloth

Created by Brant Faircloth on 22 October 2011 14:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: 

"""
import os
import sys
import argparse

# from estimate_p_i.py, possible TODO: create common library
class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def get_output_type(name):
    """get extension from filename"""
    ext = os.path.splitext(name)[1].lstrip('.').lower()
    assert ext in ['pdf','png','tiff'], "Filetype must be one of pdf, png " +\
        "or tiff"
    return ext

def get_ints_from_items(string, name = 'time'):
    """Convert times input as string to a list"""
    try:
        times = [int(i) for i in string.split(',')]
    except Exception:
        raise argparse.ArgumentTypeError("Cannot convert {} to list of " +\
                "integers".format(name))
    return times

def get_strings_from_items(string, name = 'locus'):
    """Convert times input as string to a list"""
    try:
        times = [str(i) for i in string.split(',')]
    except Exception:
        raise argparse.ArgumentTypeError("Cannot convert {} to list of " +\
                "loci".format(name))
    return times

def get_list_from_ranges(string):
    """Convert ranges entered as string to nested list"""
    try:
        ranges = [[int(j) for j in i.split('-')] for i in string.split(',')]
    except:
        raise argparse.ArgumentTypeError("Cannot convert spans to list of integers")
    return ranges
