"""
Module defining basal coancestry matrix factory interfaces and associated error checking routines.
"""

from typing import Any
from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory

from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix


class DenseCoancestryMatrixFactory(CoancestryMatrixFactory):
    """
    Factory class for producing CoancestryMatrix objects.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseCoancestryMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseCoancestryMatrixFactory, self).__init__(**kwargs)

    ############################## Object Methods ##############################
    def from_gmat(
            self, 
            gmat: GenotypeMatrix, 
            **kwargs: dict
        ) -> DenseCoancestryMatrix:
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
        out : DenseCoancestryMatrix
            A dense coancestry matrix.
        """
        return DenseCoancestryMatrix.from_gmat(
            gmat = gmat,
            **kwargs
        )




################################################################################
################################## Utilities ###################################
################################################################################
def check_is_DenseCoancestryMatrixFactory(v: object, vname: str) -> None:
    """
    Check if object is of type DenseCoancestryMatrixFactory. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseCoancestryMatrixFactory):
        raise TypeError("variable '{0}' must be a DenseCoancestryMatrixFactory".format(vname))
