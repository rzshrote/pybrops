"""
Module defining interfaces and associated error checking routines for matrices
that can have variant positions placed on a genetic map.
"""

__all__ = [
    "GeneticMappableMatrix",
    "check_is_GeneticMappableMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from pybrops.core.mat.VariantMatrix import VariantMatrix
from pybrops.popgen.gmap.GeneticMap import GeneticMap
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction

class GeneticMappableMatrix(
        VariantMatrix,
        metaclass = ABCMeta,
    ):
    """
    Abstract class for variant matrices that can be interpolated using a
    GeneticMap.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################

    ################# Interpolation Methods ################
    @abstractmethod
    def interp_genpos(
            self, 
            gmap: GeneticMap, 
            **kwargs: dict
        ) -> None:
        """
        Interpolate genetic map postions for variants using a GeneticMap

        Parameters
        ----------
        gmap : GeneticMap
            A genetic map from which to interopolate genetic map postions for
            loci within the VariantMatrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def interp_xoprob(
            self, 
            gmap: GeneticMap, 
            gmapfn: GeneticMapFunction, 
            **kwargs: dict
        ) -> None:
        """
        Interpolate genetic map positions AND crossover probabilities between
        sequential markers using a GeneticMap and a GeneticMapFunction.

        Parameters
        ----------
        gmap : GeneticMap
            A genetic map from which to interopolate genetic map postions for
            loci within the VariantMatrix.
        gmapfn : GeneticMapFunction
            A genetic map function from which to interpolate crossover
            probabilities for loci within the VariantMatrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_GeneticMappableMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type GeneticMappableMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GeneticMappableMatrix):
        raise TypeError("'{0}' must be a GeneticMappableMatrix".format(vname))
