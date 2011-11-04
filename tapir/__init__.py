"""
File: __init__.py
Author: Brant Faircloth

Created by Brant Faircloth on 22 October 2011 14:10 PDT (-0700)
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

"""

__author__ = 'Brant C. Faircloth, Jonathan Chang, Mike E. Alfaro'
__copyright__ = 'Copyright (c) 2011, Brant C. Faircloth, Jonathan Chang, Mike E. Alfaro'
__credits__ = ['Brant Faircloth', 'Jonathan Chang','Mike Alfaro']
__license__ = 'http://www.opensource.org/licenses/BSD-3-Clause'
__version__ = '1.0'
__maintainer__ = 'Brant Faircloth'
__email__ = 'brant dot faircloth at gmail dot com'
__status__ = 'testing'

import sys

if sys.version_info < (2, 7):
    raise "must use Python 2.7 or greater" # pragma: no cover

try: # pragma: no cover
    import numpy
    assert numpy.version.version > '1.3', "tapir requires >= numpy 1.3"
except ImportError: # pragma: no cover
    raise ImportError('numpy does not seem to be installed. Please see the user guide.') # pragma: no cover

try: # pragma: no cover
    import scipy
    assert scipy.version.version > '0.9', "tapir requires >= scipy 0.9.0"
except ImportError: # pragma: no cover
    raise ImportError('scipy does not seem to be installed. Please see the user guide.') # pragma: no cover

try:
    import dendropy
except ImportError: # pragma: no cover
    raise ImportError('dendropy does not seem to be installed. Please see the user guide.') # pragma: no cover

from tests import test

from db import *
from base import *
from compute import *
from pkg_resources import resource_filename

def get_hyphy_conf():
    return resource_filename(__name__, 'data/models_and_rates.bf')

def get_test_files():
    return resource_filename(__name__, 'tests/test-data')
