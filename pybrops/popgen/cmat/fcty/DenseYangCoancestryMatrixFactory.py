"""
Module defining basal coancestry matrix factory interfaces and associated error checking routines.
"""

from typing import Any, Union

import numpy
from pybrops.popgen.cmat.DenseYangCoancestryMatrix import DenseYangCoancestryMatrix
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory

from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix


class DenseYangCoancestryMatrixFactory(CoancestryMatrixFactory):
    """
    Factory class for producing CoancestryMatrix objects.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseYangCoancestryMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseYangCoancestryMatrixFactory, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def from_gmat(
            self, 
            gmat: GenotypeMatrix, 
            p_anc: Union[numpy.ndarray,float,None] = None, 
            **kwargs: dict
        ) -> DenseYangCoancestryMatrix:
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
        out : DenseYangCoancestryMatrix
            A dense coancestry matrix.
        """
        return DenseYangCoancestryMatrix.from_gmat(
            gmat = gmat,
            p_anc = p_anc,
            **kwargs
        )




################################################################################
################################## Utilities ###################################
################################################################################
def check_is_DenseYangCoancestryMatrixFactory(v: Any, vname: str) -> None:
    """
    Check if object is of type DenseYangCoancestryMatrixFactory. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseYangCoancestryMatrixFactory):
        raise TypeError("variable '{0}' must be a DenseYangCoancestryMatrixFactory".format(vname))
