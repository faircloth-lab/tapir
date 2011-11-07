Purpose
*******

*tapir* contains programs to estimate and plot phylogenetic informativeness for
large datasets.


Citing tapir
************

When using tapir, please cite:

- Faircloth BC, Chang J, Alfaro ME: *tapir* enables high throughput analysis of
  phylogenetic informativeness. `<https://github.com/faircloth-lab/tapir>`_

- Townsend JP: Profiling phylogenetic informativeness. Systematic Biol. 2007,
  56:222-231.

- Pond SLK, Frost SDW, Muse SV: HyPhy: hypothesis testing using phylogenies.
  Bioinformatics 2005, 21:676-679.

Dependencies
************

- hyphy2 (please `download <https://github.com/faircloth-lab/tapir/downloads>`_ or build a **single-threaded** hyphy2)
- Python 2.6
- numpy
- scipy
- dendropy

Installation
************
At the moment, the easiest way to install the program is::

    git clone git://github.com/faircloth-lab/tapir.git /path/to/tapir

To run tests::

    cd /path/to/tapir/
    python test/test_townsend_code.py

Use
***

The `estimate_p_i.py` code calls a batch file for hyphy that is in
`templates/`.  This file needs to be in the same position relative to
wherever you put `estimate_p_i.py`.  If you install thins as above, you'll
be fine, for the moment.

To run::

    cd /path/to/tapir/

    python tapir_compute.py Input_Folder_of_Nexus_Files/ Input.tree \
        --output Output_Directory \
        --epochs=32-42,88-98,95-105,164-174 \
        --times=37,93,100,170 \
        --multiprocessing

`--multiprocessing` is optional, without it, each locus will be run
consecutively.

If you have already run the above and saved results to your output
folder (see below), you can use the pre-existing site-rate records
rather than estimating those again with::

     python tapir_compute.py Input_Folder_of_Site_Rate_JSON_Files/ Input.tree \
        --output Output_Directory \
        --epochs=32-42,88-98,95-105,164-174 \
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
`wiki <https://github.com/faircloth-lab/tapir/wiki/getting-data-from-the-database>`_ 
for examples.

Acknowledgements
****************

BCF thanks SP Hubbell, PA Gowaty, RT Brumfield, TC Glenn, NG Crawford,
JE McCormack, and M Reasel. JHLC and MEA thank J Eastman and J Brown for
thoughtful comments about PI. We thank Francesc Lopez-Giraldez and
Jeffrey Townsend for providing us with a copy of their web-application
source code and helpful discussion.
