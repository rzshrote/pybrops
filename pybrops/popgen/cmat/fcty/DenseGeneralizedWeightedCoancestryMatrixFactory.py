"""
Module defining basal coancestry matrix factory interfaces and associated error checking routines.
"""

from numbers import Real
from typing import Any, Union
import numpy
from pybrops.popgen.cmat.DenseGeneralizedWeightedCoancestryMatrix import DenseGeneralizedWeightedCoancestryMatrix
from pybrops.popgen.cmat.fcty.DenseCoancestryMatrixFactory import DenseCoancestryMatrixFactory
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix


class DenseGeneralizedWeightedCoancestryMatrixFactory(DenseCoancestryMatrixFactory):
    """
    Factory class for producing CoancestryMatrix objects.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseGeneralizedWeightedCoancestryMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseGeneralizedWeightedCoancestryMatrixFactory, self).__init__(**kwargs)

    ############################## Object Methods ##############################
    def from_gmat(
            self, 
            gmat: GenotypeMatrix, 
            mkrwt: Union[numpy.ndarray,Real,None] = None,
            afreq: Union[numpy.ndarray,Real,None] = None, 
            **kwargs: dict
        ) -> DenseGeneralizedWeightedCoancestryMatrix:
        """
        Create a CoancestryMatrix from a GenotypeMatrix.

        Parameters
        ----------
        gmat : GenotypeMatrix
            Input genotype matrix from which to calculate coancestry.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseGeneralizedWeightedCoancestryMatrix
            A dense coancestry matrix.
        """
        return DenseGeneralizedWeightedCoancestryMatrix.from_gmat(
            gmat = gmat,
            mkrwt = mkrwt,
            afreq = afreq,
            **kwargs
        )




################################################################################
################################## Utilities ###################################
################################################################################
def check_is_DenseGeneralizedWeightedCoancestryMatrixFactory(v: object, vname: str) -> None:
    """
    Check if object is of type DenseGeneralizedWeightedCoancestryMatrixFactory. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseGeneralizedWeightedCoancestryMatrixFactory):
        raise TypeError("variable '{0}' must be a DenseGeneralizedWeightedCoancestryMatrixFactory".format(vname))
