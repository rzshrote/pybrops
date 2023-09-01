"""
Module implementing the extraction of true breeding value.
"""

from typing import Union
import numpy
from pybrops.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import check_is_PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

class TrueBreedingValue(BreedingValueProtocol):
    """
    Class implementing the extraction of true breeding value.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            gpmod: GenomicModel, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class TrueBreedingValue.

        Parameters
        ----------
        gpmod : GenomicModel
            A true genomic prediction model from which to calculate true
            breeding values.
        """
        super(TrueBreedingValue, self).__init__(**kwargs)
        self.gpmod = gpmod

    ############################ Object Properties #############################

    ############### Genomic Model Properties ###############
    @property
    def gpmod(self) -> GenomicModel:
        """Genomic prediction model."""
        return self._gpmod
    @gpmod.setter
    def gpmod(self, value: GenomicModel) -> None:
        """Set genomic prediction model."""
        check_is_GenomicModel(value, "gpmod")
        self._gpmod = value

    ############################## Object Methods ##############################
    def estimate(
            self, 
            ptobj: Union[PhenotypeDataFrame,BreedingValueMatrix,numpy.ndarray], 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            miscout: dict = None, 
            **kwargs: dict
        ) -> BreedingValueMatrix:
        """
        Estimate breeding values.

        Parameters
        ----------
        ptobj : BreedingValueMatrix, PhenotypeDataFrame, numpy.ndarray
            An object containing phenotype data. Must be a matrix of breeding
            values or a phenotype data frame.
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        gpmod : GenomicModel
            Genomic model used for predicting genotypes.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : BreedingValueMatrix
            A matrix of breeding values.
        """
        # calculate true breeding values
        bvmat = self.gpmod.gebv(gtobj)

        return bvmat
