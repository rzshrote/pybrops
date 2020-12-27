class GeneticMapFunction:
    """docstring for GeneticMapFunction."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(GeneticMapFunction, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ########## Mapping & Inverse Mapping Methods ###########
    def mapfn(self, d):
        """
        Convert genetic distance (d) to recombination probability (r).

        Parameters
        ----------
        d : numpy.ndarray
            An array of genetic map distances in Morgans.

        Returns
        -------
        r : numpy.ndarray
            An array of recombination probabilities for each corresponding
            distance.
        """
        raise NotImplementedError("method is abstract")

    def invmapfn(self, r):
        """
        convert recombination probability (r) to genetic distance (d).

        Parameters
        ----------
        r : numpy.ndarray
            An array of recombination probabilities.

        Returns
        -------
        d : numpy.ndarray
            An array of genetic distances in Morgans.
        """
        raise NotImplementedError("method is abstract")

    ########## Recombination Probability Methods ###########
    def rprob1g(self, gmap, vrnt_chrgrp, vrnt_genpos):
        """
        Calculate sequential recombination probabilities using genetic distances.
        Calculate recombination probabilities between successive entries along
        a chromosome. Supply 0.5 across chromosomes.
        Example:
            vrnt_chrgrp = [1,  1,  2,  2,  2,  3,  3]
            vrnt_genpos = ...
            xo         = [.5, .2, .2, .1, .3, .5, .2]

        Parameters
        ----------
        gmap : GeneticMap
            GeneticMap object for calculating genetic distances between successive
            entries along a chromosome.
        vrnt_chrgrp : numpy.ndarray
            An array assigning chromosomes to groups. Must be sorted.
        vrnt_genpos : numpy.ndarray
            An array of genetic positions. Must be sorted and correspond with
            vrnt_chrgrp.

        Returns
        -------
        r : numpy.ndarray
        """
        raise NotImplementedError("method is abstract")

    def rprob2g(self, gmap, vrnt_chrgrp, vrnt_genpos):
        """
        Calculate pairwise recombination probabilities using genetic distances.
        Calculate a recombination probability matrix.

        Parameters
        ----------
        vrnt_chrgrp : numpy.ndarray
            An array assigning chromosomes to groups. Must be sorted.
        vrnt_genpos : numpy.ndarray
            An array of genetic positions. Must be sorted and correspond with
            vrnt_chrgrp.
        key : tuple
            A tuple of array calculation regions. Supply these if only a section
            of the final recombination matrix is desired.

        Returns
        -------
        r : numpy.ndarray
            A 2D array of recombination probabilities.
        """
        raise NotImplementedError("method is abstract")

    def rprob1p(self, gmap, vrnt_chrgrp, vrnt_phypos):
        """
        Calculate sequential recombination probabilities using physical distances.
        """
        raise NotImplementedError("method is abstract")

    def rprob2p(self, gmap, vrnt_chrgrp, vrnt_phypos):
        """
        Calculate pairwise recombination probabilities using physical distances.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GeneticMapFunction(v):
    return isinstance(v, GeneticMapFunction)

def check_is_GeneticMapFunction(v, varname):
    if not isinstance(v, GeneticMapFunction):
        raise TypeError("'%s' must be an GeneticMapFunction." % varname)

def cond_check_is_GeneticMapFunction(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GeneticMapFunction(v, varname)
