"""
Module defining basal coancestry matrix factory interfaces and associated error checking routines.
"""

__all__ = [
    "CoancestryMatrixFactory",
    "check_is_CoancestryMatrixFactory"
]

from pybrops.popgen.cmat.CoancestryMatrix import CoancestryMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class CoancestryMatrixFactory:
    """
    Factory class for producing CoancestryMatrix objects.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for CoancestryMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(CoancestryMatrixFactory, self).__init__(**kwargs)

    ############################## Object Methods ##############################
    def from_gmat(
            cls, 
            gmat: GenotypeMatrix, 
            **kwargs: dict
        ) -> CoancestryMatrix:
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
        out : CoancestryMatrix
            A coancestry matrix.
        """
        raise NotImplementedError("class method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_CoancestryMatrixFactory(v: object, vname: str) -> None:
    """
    Check if object is of type CoancestryMatrixFactory. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, CoancestryMatrixFactory):
        raise TypeError("variable '{0}' must be a CoancestryMatrixFactory".format(vname))
