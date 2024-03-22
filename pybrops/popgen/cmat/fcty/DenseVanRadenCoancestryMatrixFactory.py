"""
Module defining basal coancestry matrix factory interfaces and associated error checking routines.
"""

__all__ = [
    "DenseVanRadenCoancestryMatrixFactory",
    "check_is_DenseVanRadenCoancestryMatrixFactory",
]

from typing import Union

import numpy
from pybrops.popgen.cmat.DenseVanRadenCoancestryMatrix import DenseVanRadenCoancestryMatrix
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class DenseVanRadenCoancestryMatrixFactory(
        CoancestryMatrixFactory,
    ):
    """
    Factory class for producing CoancestryMatrix objects.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseVanRadenCoancestryMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseVanRadenCoancestryMatrixFactory, self).__init__(**kwargs)

    ############################## Object Methods ##############################
    def from_gmat(
            self, 
            gmat: GenotypeMatrix, 
            p_anc: Union[numpy.ndarray,float,None] = None, 
            **kwargs: dict
        ) -> DenseVanRadenCoancestryMatrix:
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
        out : DenseVanRadenCoancestryMatrix
            A dense coancestry matrix.
        """
        return DenseVanRadenCoancestryMatrix.from_gmat(
            gmat = gmat,
            p_anc = p_anc,
            **kwargs
        )



################################## Utilities ###################################
def check_is_DenseVanRadenCoancestryMatrixFactory(v: object, vname: str) -> None:
    """
    Check if object is of type DenseVanRadenCoancestryMatrixFactory. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseVanRadenCoancestryMatrixFactory):
        raise TypeError("variable '{0}' must be a DenseVanRadenCoancestryMatrixFactory".format(vname))
