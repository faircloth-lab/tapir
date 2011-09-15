    ***************************************************
    *                                                 *
    * pd-ht:  PhyDesign - High-Throughput             *
    *                                                 *
    * (c) 2011 Brant Faircloth, Jonathan Chang,       *
    * Mike Alfaro                                     *
    *                                                 *
    * PhyDesign was created by the Townsend Lab       *
    * (http://phydesign.townsend.yale.edu)            *
    *                                                 *
    * To cite Phydesign, please use:                  *
    *                                                 *
    *   - J.P. Townsend, 20007. Profiling             *
    *     phylogenetic informativeness. Systematic    *
    *     Biology, 56(2), 222-231.                    *
    *                                                 *
    *   - Pond, S.L.K., Frost, S.D.W., and S.V. Muse, *
    *     2005. Hyphy: hypothesis testing using       *
    *     phylogenies. Bioinformatics, 21(5), 676-9.  *
    *                                                 *
    * Many thanks to Francesc Lopez-Giraldez and      *
    * Jeffrey Townsend for providing us with a copy   *
    * of their web-application source code.           *
    *                                                 *
    ***************************************************

## Dependencies

 * hyphy (please [download](https://github.com/downloads/BadDNA/pd-ht/hyphy2.osx.gz) or build the single-threaded hyphy2)
 * Python 2.7
 * numpy
 * scipy
 * dendropy

## Installation

 * git clone git://github.com/BadDNA/pd-ht.git /path/to/pd-ht

## Use

The `estimate_p_i.py` code calls a batch file for hyphy that is in
`templates/`.  This file needs to be in the same position relative to
wherever you put `estimate_p_i.py`.

    cd /path/to/pd-ht/

    python estimate_p_i.py Input_Folder_of_Nexus_Files/ Input.tree \
        --output Output_Directory \
        --epochs=32-42,88-98,95-105,164-174 \
        --times=37,93,100,170 \
        --multiprocessing

`--multiprocessing` is optional, without it, each locus will be run
consecutively.

If you have already run the above and saved results to your output
folder (see below), you can use the pre-existing site-rate records
rather than estimating those again with:

     python estimate_p_i.py Input_Folder_of_Site_Rate_JSON_Files/ Input.tree \
        --output Output_Directory \
        --epochs=32-42,88-98,95-105,164-174 \
        --times=37,93,100,170 \
        --multiprocessing \
        --site-rates

## Results

Results will be in an [sqlite](http://www.sqlite.org/) database in the
output directory of your choosing.  This directory also holds site rate
files in [JSON](http://www.json.org/) format for each locus passed
through `estimate_p_i.py`.
