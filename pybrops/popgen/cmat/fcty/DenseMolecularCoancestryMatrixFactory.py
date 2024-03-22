"""
Module defining basal coancestry matrix factory interfaces and associated error checking routines.
"""

__all__ = [
    "DenseMolecularCoancestryMatrixFactory",
    "check_is_DenseMolecularCoancestryMatrixFactory",
]

from pybrops.popgen.cmat.DenseMolecularCoancestryMatrix import DenseMolecularCoancestryMatrix
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class DenseMolecularCoancestryMatrixFactory(
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
        Constructor for DenseMolecularCoancestryMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseMolecularCoancestryMatrixFactory, self).__init__(**kwargs)

    ############################## Object Methods ##############################
    def from_gmat(
            self, 
            gmat: GenotypeMatrix, 
            **kwargs: dict
        ) -> DenseMolecularCoancestryMatrix:
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
        out : DenseMolecularCoancestryMatrix
            A dense coancestry matrix.
        """
        return DenseMolecularCoancestryMatrix.from_gmat(
            gmat = gmat,
            **kwargs
        )



################################## Utilities ###################################
def check_is_DenseMolecularCoancestryMatrixFactory(v: object, vname: str) -> None:
    """
    Check if object is of type DenseMolecularCoancestryMatrixFactory. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseMolecularCoancestryMatrixFactory):
        raise TypeError("variable '{0}' must be a DenseMolecularCoancestryMatrixFactory".format(vname))
