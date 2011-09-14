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

## Use

    python estimate_p_i.py Input_Folder_of_Nexus_files/ Input.tree \
        --output Output_Directory \
        --epochs=32-42,88-98,95-105,164-174 \
        --times=37,93,100,170 \
        --multiprocessing
