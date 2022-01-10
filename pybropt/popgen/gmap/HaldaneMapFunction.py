import numpy

from pybropt.popgen.gmap.GeneticMapFunction import GeneticMapFunction

class HaldaneMapFunction(GeneticMapFunction):
    """docstring for HaldaneMapFunction."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Create a Haldane mapping function object.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(HaldaneMapFunction, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ########## Mapping & Inverse Mapping Methods ###########
    def mapfn(self, d):
        """
        Convert genetic distances in Morgans to recombination probabilities using
        the Haldane mapping function (Haldane, 1919).

        This is a bare-bones function and does no parameter checking for errors.

        Parameters
        ----------
        d : numpy.ndarray
            An array of genetic distances in Morgans. This can be an array of any
            shape.

        Returns
        -------
        r : numpy.ndarray
            An array of recombination probabilities. The shape of the array is
            the same shape as that of 'd'.
        """
        # convert d to r
        r = 0.5 * (1.0 - numpy.exp(-2.0 * d))
        return r

    def invmapfn(self, r):
        """
        Convert recombination probabilities between two loci to genetic map
        distances using the Haldane mapping function (Haldane, 1919).

        This is a bare-bones function and does no parameter checking for errors.

        Parameters
        ----------
        r : numpy.ndarray
            An array of recombination probabilities between two loci.

        Returns
        -------
        d : numpy.ndarray
            An array of genetic map distances as defined by the Haldane mapping
            function.
        """
        # convert r to d
        d = -0.5 * numpy.log(1.0 - (2.0 * r))
        return d

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
        return self.mapfn(gmap.gdist1g(vrnt_chrgrp, vrnt_genpos))

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
        return self.mapfn(gmap.gdist2g(vrnt_chrgrp, vrnt_genpos))

    def rprob1p(self, gmap, vrnt_chrgrp, vrnt_phypos):
        """
        Calculate sequential recombination probabilities using physical distances.
        """
        return self.mapfn(gmap.gdist1p(vrnt_chrgrp, vrnt_phypos))

    def rprob2p(self, gmap, vrnt_chrgrp, vrnt_phypos):
        """
        Calculate pairwise recombination probabilities using physical distances.
        """
        return self.mapfn(gmap.gdist2p(vrnt_chrgrp, vrnt_phypos))



################################################################################
################################## Utilities ###################################
################################################################################
def is_HaldaneMapFunction(v):
    """Return whether an object is a HaldaneMapFunction or not"""
    return isinstance(v, HaldaneMapFunction)

def check_is_HaldaneMapFunction(v, vname):
    """Raise TypeError if object is not a HaldaneMapFunction"""
    if not isinstance(v, HaldaneMapFunction):
        raise TypeError("variable '{0}' must be a HaldaneMapFunction".format(vname))

def cond_check_is_HaldaneMapFunction(v, vname, cond=(lambda s: s is not None)):
    """If object is not None, raise TypeError if object is not a HaldaneMapFunction"""
    if cond(v):
        check_is_HaldaneMapFunction(v, vname)
