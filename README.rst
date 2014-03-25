Purpose
*******

*tapir* contains programs to estimate and plot phylogenetic informativeness for
large datasets.


Citing tapir
************

When using tapir, please cite:

- Faircloth BC, Chang J, Alfaro ME: *TAPIR* enables high throughput analysis of
  phylogenetic informativeness. `arXiv:1202.1215 <http://arxiv.org/abs/1202.1215>`_

- Townsend JP: Profiling phylogenetic informativeness. Systematic Biol. 2007,
  56:222-231.

- Pond SLK, Frost SDW, Muse SV: HyPhy: hypothesis testing using phylogenies.
  Bioinformatics 2005, 21:676-679.

Dependencies
************

- hyphy2 (please download [`osx <http://s3.faircloth-lab.org/packages/hyphy2.osx.gz>`_ or `linux <http://s3.faircloth-lab.org/packages/hyphy2.linux.gz>`_] or build a **single-threaded** hyphy2)
- Python 2.6
- numpy
- scipy
- dendropy

Installation
*************

For **ALL** platforms, you must download a hyphy binary for your platform (osx
or linux) and place that within your $PATH.  To install the other dependencies
(numpy, scipy), you may need to install a Fortran compiler on linux/osx, which
is a dependency of numpy and scipy.

Linux
------

Install hyphy::

    wget https://github.com/downloads/faircloth-lab/tapir/hyphy2.linux.gz
    gunzip hyphy2.*.gz
    chmod 0700 hyphy2.*
    mv hyphy2.* ~/Bin/hyphy2

On linux (ubuntu/debian), use::

    apt-get install gfortran libatlas-base-dev liblapack-dev

Install tapir and dependencies, which include numpy and scipy (the
reason we installed the dependencies above)::

    pip install tapir

To plot results, you will also need to::

    apt-get install r-base r-base-dev
    pip install rpy2

OSX
---

Install hyphy::

    wget https://github.com/downloads/faircloth-lab/tapir/hyphy2.osx.gz
    gunzip hyphy2.*.gz
    chmod 0700 hyphy2.*
    mv hyphy2.* ~/Bin/hyphy2

It is easiest just to install the `scipy superpack
<http://stronginference.com/scipy-superpack/>`_.  This will install the
dependencies that tapir needs.  After installing the superpack, using
pip, install tapir::

    pip install tapir

Alternatively, you can simply try to install tapir using::

    pip install tapir

To plot results, you need to `install R
<http://cran.r-project.org/bin/macosx/>`_ and then install rpy2::

    pip install rpy2


Other OSs
----------

Download the source for hyphy from `<http://www.datam0nk3y.org/hyphy/doku.php?id=start>`_
and compile a **single-threaded** version.

Install `numpy <http://numpy.scipy.org>`_, `scipy <http://scipy.org>`_,
and `dendropy <http://packages.python.org/DendroPy/>`_ for your
platform.  Alternatively, you can try::

    pip install tapir

Which may install the necessary dependencies.  Alternatively, install the 
dependenceis manually, and then::

    wget http://pypi.python.org/packages/source/t/tapir/tapir-1.0.tar.gz
    tar -xzvf tapir-1.0.tar.gz
    cd tapir*
    python setup.py build
    python setup.py test
    python setup.py install


Testing
*******

If you didn't run the tests using `python setup.py test` above, you can also::

    import tapir
    tapir.test()

Use
***
To run::

    python tapir_compute.py Input_Folder_of_Nexus_Files/ Input.tree \
        --output Output_Directory \
        --intervals=32-42,88-98,95-105,164-174 \
        --times=37,93,100,170 \
        --multiprocessing

`--multiprocessing` is optional, without it, each locus will be run
consecutively.

If you have already run the above and saved results to your output
folder (see below), you can use the pre-existing site-rate records
rather than estimating those again with::

     python tapir_compute.py Input_Folder_of_Site_Rate_JSON_Files/ Input.tree \
        --output Output_Directory \
        --intervals=32-42,88-98,95-105,164-174 \
        --times=37,93,100,170 \
        --multiprocessing \
        --site-rates

Results
*******

tapir writes results to a `sqlite <http://www.sqlite.org/>`_ database in the
output directory of your choosing.  This directory also holds site rate
files in `JSON <http://www.json.org/>`_ format for each locus passed
through `tapir_compute.py`.

You can access the results in the database as follows.  For more examples,
including plotting, see the 
`documentation <http://faircloth-lab.github.com/tapir/>`_

- crank up sqlite::

    sqlite3  Output_Directory/phylogenetic-informativeness.sqlite

- get integral data for all epochs::

    select locus, interval, pi from loci, interval where loci.id = interval.id

- get integral data for a specific epoch::

    select locus, interval, pi from loci, interval 
    where interval = '95-105' and loci.id = interval.id;

- get the count of loci having max(PI) at different epochs::

    create temporary table max as select id, max(pi) as max from interval group by id;

    create temporary table t as select interval.id, interval, max from interval, max 
    where interval.pi = max.max;

    select interval, count(*) from t group by interval;

Plotting Results
****************

tapir contains plotting scripts to help you plot data within a results database
and compare data between different databases.  tapir uses RPY and R to
do this.  You can also plot data directly in R.  Until we finish the
documentation, please see the 
`wiki <https://github.com/faircloth-lab/tapir/wiki/getting-data-from-the-database(s)>`_ 
for examples.

Acknowledgements
****************

BCF thanks SP Hubbell, PA Gowaty, RT Brumfield, TC Glenn, NG Crawford,
JE McCormack, and M Reasel. JHLC and MEA thank J Eastman and J Brown for
thoughtful comments about PI. We thank Francesc Lopez-Giraldez and
Jeffrey Townsend for providing us with a copy of their web-application
source code and helpful discussion.
