"""
Module implementing the Haldane genetic map function and associated error checking routines.
"""

from typing import Any
import numpy
from pybrops.popgen.gmap.GeneticMap import GeneticMap
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction

class HaldaneMapFunction(GeneticMapFunction):
    """
    A concrete class for the Haldane genetic map function.

    The purpose of this concrete class is to implement functionality for:
        1) Converting genetic distance to recombination probability.
        2) Converting recombination probability to genetic distance.
        3) Converting physical distance to recombination probability.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for a Haldane mapping function object.

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
    def mapfn(
            self, 
            d: numpy.ndarray
        ) -> numpy.ndarray:
        """
        Convert genetic distances in Morgans to recombination probabilities using
        the Haldane mapping function (Haldane, 1919).

        This is a bare-bones function and does no parameter checking for errors.

        Parameters
        ----------
        d : numpy.ndarray
            An array of genetic distances in Morgans. This can be an array of
            any shape.

        Returns
        -------
        r : numpy.ndarray
            An array of recombination probabilities. The shape of the array is
            the same shape as that of ``d``.
        """
        # convert d to r
        r = 0.5 * (1.0 - numpy.exp(-2.0 * d))
        return r

    def invmapfn(
            self, 
            r: numpy.ndarray
        ) -> numpy.ndarray:
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
            function. The shape of the array is the same shape as that of ``r``.
        """
        # convert r to d
        d = -0.5 * numpy.log(1.0 - (2.0 * r))
        return d

    ########## Recombination Probability Methods ###########
    def rprob1g(
            self, 
            gmap: GeneticMap, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray
        ) -> numpy.ndarray:
        """
        Calculate sequential recombination probabilities using genetic distances.
        Calculate recombination probabilities between successive entries along
        a chromosome. Supply 0.5 across chromosomes.

        Example::

            vrnt_chrgrp = [ 1,  1,  2,  2,  2,  3,  3]
            vrnt_genpos = [...]
            xo          = [.5, .2, .5, .1, .3, .5, .2]

        Parameters
        ----------
        gmap : GeneticMap
            GeneticMap object for calculating genetic distances between successive
            entries along a chromosome.
        vrnt_chrgrp : numpy.ndarray
            An array assigning chromosomes to groups. Must be sorted.
        vrnt_genpos : numpy.ndarray
            An array of genetic positions. Must be sorted and correspond with
            ``vrnt_chrgrp``.

        Returns
        -------
        r : numpy.ndarray
            A 1D array of recombination probabilities.
        """
        return self.mapfn(gmap.gdist1g(vrnt_chrgrp, vrnt_genpos))

    def rprob2g(
            self, 
            gmap: GeneticMap, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray
        ) -> numpy.ndarray:
        """
        Calculate pairwise recombination probabilities using genetic distances.
        Calculate a recombination probability matrix.

        Parameters
        ----------
        gmap : GeneticMap
            GeneticMap object for calculating genetic distances between successive
            entries along a chromosome.
        vrnt_chrgrp : numpy.ndarray
            An array assigning chromosomes to groups. Must be sorted.
        vrnt_genpos : numpy.ndarray
            An array of genetic positions. Must be sorted and correspond with
            ``vrnt_chrgrp``.

        Returns
        -------
        r : numpy.ndarray
            A 2D array of recombination probabilities.
        """
        return self.mapfn(gmap.gdist2g(vrnt_chrgrp, vrnt_genpos))

    def rprob1p(
            self, 
            gmap: GeneticMap, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray
        ) -> numpy.ndarray:
        """
        Calculate sequential recombination probabilities using physical distances.

        Parameters
        ----------
        gmap : GeneticMap
            GeneticMap object for calculating genetic distances between successive
            entries along a chromosome.
        vrnt_chrgrp : numpy.ndarray
            An array assigning chromosomes to groups. Must be sorted.
        vrnt_phypos : numpy.ndarray
            An array of physical positions. Must be sorted and correspond with
            ``vrnt_chrgrp``.

        Returns
        -------
        r : numpy.ndarray
            A 1D array of recombination probabilities.
        """
        return self.mapfn(gmap.gdist1p(vrnt_chrgrp, vrnt_phypos))

    def rprob2p(
            self, 
            gmap: GeneticMap, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray
        ) -> numpy.ndarray:
        """
        Calculate pairwise recombination probabilities using physical distances.

        Parameters
        ----------
        gmap : GeneticMap
            GeneticMap object for calculating genetic distances between successive
            entries along a chromosome.
        vrnt_chrgrp : numpy.ndarray
            An array assigning chromosomes to groups. Must be sorted.
        vrnt_phypos : numpy.ndarray
            An array of physical positions. Must be sorted and correspond with
            ``vrnt_chrgrp``.

        Returns
        -------
        r : numpy.ndarray
            A 2D array of recombination probabilities.
        """
        return self.mapfn(gmap.gdist2p(vrnt_chrgrp, vrnt_phypos))



################################################################################
################################## Utilities ###################################
################################################################################
def is_HaldaneMapFunction(v: Any) -> bool:
    """
    Determine whether an object is a HaldaneMapFunction.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a HaldaneMapFunction object instance.
    """
    return isinstance(v, HaldaneMapFunction)

def check_is_HaldaneMapFunction(v: Any, varname: str) -> None:
    """
    Check if object is of type HaldaneMapFunction. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, HaldaneMapFunction):
        raise TypeError("'{0}' must be of type HaldaneMapFunction.".format(varname))
