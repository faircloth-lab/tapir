.. _installation:

**************
Installation
**************

picme_ has several dependencies.  At a minimum you will need to install:

- `hyphy 2 <http://www.datam0nk3y.org/hyphy/>`_ (must be
  single-threaded)
- `numpy 1.3.x <http://numpy.scipy.org>`_
- `scipy 0.9.0 <http://scipy.org>`_
- `dendropy 3.9.0 <http://packages.python.org/DendroPy/>`_

Installing hyphy
================

picme_ uses hyphy_ for substitution model selection and site rate
estimation.  Additionally, hyphy needs to be compiled as a
single-threaded binary.  We provide appropriate binaries for:

- `os x 10.7 <https://github.com/downloads/faircloth-lab/picme/hyphy2.osx.gz>`_
- linux

You need to download these and place them within your `$PATH`.  If this
statement confuses you, then you should `read this article on UNIX paths
<http://kb.iu.edu/data/acar.html>`_.

Installing numpy
================

You first need to install numpy_.  Barring platform-specific options, you can 
accomplish this by running::

    easy_install numpy

Sometimes on OSX, this is problematic.  There is a `binary installer
<http://sourceforge.net/projects/numpy/files/NumPy/1.6.1/numpy-1.6.1-py2.6-python.org-macosx10.3.dmg/download>`_
that you can use in this case.

Installing scipy
================

You also need to install scipy_.   Barring platform-specific options, you can 
accomplish this by running::

    easy_install scipy

Sometimes on OSX, this is problematic.  There is a `binary installer
<http://sourceforge.net/projects/scipy/files/scipy/0.9.0/scipy-0.9.0-py2.7-python.org-macosx10.6.dmg/download>`_
that you can use in this case.

Installing dendropy
===================

Installing dendropy_ is easy using easy_install::

    easy_install dendropy



rpy2 [optional]
=======================

To produce pretty and relatively pain-free graphics, you will need to
install rpy2_.  This is generally easily done using::

    easy_install rpy2

Although this is currently not working on OSX Lion 
(see `here <https://bitbucket.org/lgautier/rpy2/issue/78/truncated-path-to-librblasdylib-on-os-x>`_ for
temporary workaround).

matplotlib [optional]
=============================

To produce more basic figures, you can use matplotlib which is installed
using::

    easy_install matplotlib


.. _hyphy: http://www.datam0nk3y.org/hyphy/
.. _rpy2: http://rpy.sourceforge.net/rpy2.html
.. _scipy: http://scipy.org
.. _numpy: http://numpy.scipy.org
.. _dendropy: http://packages.python.org/DendroPy
.. _picme: https://github.com/faircloth-lab/picme
