"""
Module defining basal genetic map function interfaces and associated error
checking routines.
"""

__all__ = [
    "GeneticMapFunction",
    "check_is_GeneticMapFunction",
]

from abc import ABCMeta
from abc import abstractmethod
import numpy
from pybrops.popgen.gmap.GeneticMap import GeneticMap

class GeneticMapFunction(
        metaclass = ABCMeta,
    ):
    """
    An abstract class for genetic map function objects.

    The purpose of this abstract class is to define base functionality for:
        1) Converting genetic distance to recombination probability.
        2) Converting recombination probability to genetic distance.
        3) Converting physical distance to recombination probability.
    """

    ########################## Special Object Methods ##########################
    def __repr__(
            self
        ) -> str:
        """
        Return repr(self).
        
        Returns
        -------
        out : str
            A representation of the object.
        """
        return "<{0} at {1}>".format(
            type(self).__name__,
            hex(id(self))
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Mapping & Inverse Mapping Methods ###########
    @abstractmethod
    def mapfn(
            self, 
            d: numpy.ndarray
        ) -> numpy.ndarray:
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

    @abstractmethod
    def invmapfn(
            self, 
            r: numpy.ndarray
        ) -> numpy.ndarray:
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
    @abstractmethod
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
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
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_GeneticMapFunction(v: object, vname: str) -> None:
    """
    Check if object is of type GeneticMapFunction. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GeneticMapFunction):
        raise TypeError("'{0}' must be of type GeneticMapFunction.".format(vname))
