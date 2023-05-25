"""
Module implementing dense matrices that can have their variants placed on a
genetic map and associated error checking routines.
"""

from typing import Any, Optional

import numpy
from pybrops.core.mat.DenseVariantMatrix import DenseVariantMatrix
from pybrops.popgen.gmap.GeneticMap import GeneticMap, check_is_GeneticMap
from pybrops.popgen.gmap.GeneticMappableMatrix import GeneticMappableMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction, check_is_GeneticMapFunction

class DenseGeneticMappableMatrix(DenseVariantMatrix,GeneticMappableMatrix):
    """
    Concrete class for dense variant matrices that can be interpolated using a
    GeneticMap.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            vrnt_chrgrp: Optional[numpy.ndarray] = None, 
            vrnt_phypos: Optional[numpy.ndarray] = None,
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None, 
            vrnt_xoprob: Optional[numpy.ndarray] = None,
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None, 
            vrnt_hapref: Optional[numpy.ndarray] = None,
            vrnt_mask: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DenseGeneticMappableMatrix

        Parameters
        ----------
        mat : numpy.ndarray
        vrnt_chrgrp : numpy.ndarray, None
        vrnt_phypos : numpy.ndarray, None
        vrnt_name : numpy.ndarray, None
        vrnt_genpos : numpy.ndarray, None
        vrnt_xoprob : numpy.ndarray, None
        vrnt_hapgrp : numpy.ndarray, None
        vrnt_hapalt : numpy.ndarray, None
        vrnt_hapref : numpy.ndarray, None
        vrnt_mask : numpy.ndarray, None
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseGeneticMappableMatrix, self).__init__(
            mat = mat,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask
            **kwargs
        )

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################# Interpolation Methods ################
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
        """
        # check if gmap is a GeneticMap
        check_is_GeneticMap(gmap, "gmap")

        # interpolate postions
        self.vrnt_genpos = gmap.interp_genpos(self._vrnt_chrgrp, self._vrnt_phypos)

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
        """
        # check data types
        check_is_GeneticMap(gmap, "gmap")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")

        # check if self has been sorted and grouped
        if not self.is_grouped_vrnt():
            raise RuntimeError("must be grouped first before interpolation of crossover probabilities")

        # interpolate genetic positions
        self.vrnt_genpos = gmap.interp_genpos(self._vrnt_chrgrp, self._vrnt_phypos)

        # interpolate crossover probabilities
        self.vrnt_xoprob = gmapfn.rprob1g(gmap, self._vrnt_chrgrp, self._vrnt_genpos)



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseGeneticMappableMatrix(v: object) -> bool:
    """
    Determine whether an object is a DenseGeneticMappableMatrix.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseGeneticMappableMatrix object instance.
    """
    return isinstance(v, DenseGeneticMappableMatrix)

def check_is_DenseGeneticMappableMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseGeneticMappableMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseGeneticMappableMatrix):
        raise TypeError("'{0}' must be a DenseGeneticMappableMatrix".format(vname))
