"""
File: __init__.py
Author: Brant Faircloth

Created by Brant Faircloth on 22 October 2011 14:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

Description: 

"""

from db import *
from base import *
from compute import *
#from rfunctions import *
from pkg_resources import resource_string

def get_hyphy_conf():
    return resource_filename(__name__, 'data/models_and_rates.bf')
