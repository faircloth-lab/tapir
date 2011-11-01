.. _computing:

***************************************
Computing Phylogenetic Informativeness
***************************************

After you get the dependencies installed, computing the phylogenetic
informativeness (PI) is rather easy.  The typical invocation of the
program looks like:

.. code-block:: bash

    picme_compute.py /your/input/nexus/ \
        /your/data/treefile.tree \
        --output /path/to/the/output \
        --times 1,2,3,4 \
        --intervals 1-5,5-10

This wil compute the PI from a folder of nexus files at
`/your/input/nexus/`, against the tree that is located at 
`/your/data/treefile.tree`, for the the discrete times 1,2,3,4 and
across the intervals (*aka* epochs) 1-5 and 5-10.

We describe each option (`picme_compute.py -h`) in detail, below.  We 
have separated these options into those pertinent to
:ref:`initial-computation` and those for :ref:`re-computation`.

.. _initial-computation:

PI Computation Options
**********************

Positional arguments
---------------------

    **alignments**  The folder of alignments

    **tree**  The input tree

Positional arguments are passed to the program by *position*.  For
example, to send the above to `picme_compute.py`, you would type, on the
command-line

.. code-block:: bash

    picme_compute.py alignments tree

Required arguments
------------------

We have several options that are requies, which can annoy some people.
Here, we've made these options required because it makes the
command-line invocation of the program less clumsy and more clear.

--times TIMES  Comma-separated list of start times of interest (MYA)

--epochs EPOCHS  Comma-separated list of epoch ranges of interest (MYA)


Optional arguments
------------------

These are true options in the sense that they are entirely optional -
the defaults will be chose for you if you do not pass a value along with
it's corresponding option.

-h, --help  show this help message and exit

--tree-format FORMAT  The format of the tree (`nexus` or `newick`)

--output OUTPUT  The path to the output directory

--hyphy HYPHY  The path to hyphy (if not in $PATH)

--template TEMPLATE  Path to the hypy temlate file if in non-standard
  location

--threshold THRESHOLD  Minimum number of taxa without a gap for a site
  to be considered informative

--multiprocessing  Enable parallel calculation of rates

--site-rates  Use previously calculated site rates



