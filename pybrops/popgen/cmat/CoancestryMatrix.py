"""
Module defining basal coancestry matrix interfaces and associated error checking routines.
"""

from typing import Any

import numpy
from pybrops.core.mat.SquareTaxaMatrix import SquareTaxaMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class CoancestryMatrix(SquareTaxaMatrix):
    """
    An abstract class for coancestry matrices. Coancestry matrices are square.
    Coancestry matrices are related to kinship matrices in the following manner:

    ..math:
        \\mathbf{K} = \\frac{1}{2}\\mathbf{A}

    The purpose of this abstract class is to define base functionality for:
        1) Coancestry matrix value calculation.
        2) Coancestry matrix value access.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for CoancestryMatrix class.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments for dependency injection.
        """
        super(CoancestryMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############### Genotype Data Properites ###############
    # gmat should be implemented in a sparse version of this matrix.

    ############## Coancestry Data Properites ##############
    # access using mat (inherited from Matrix)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################## Matrix conversion ###################
    def mat_asformat(self, format: str) -> numpy.ndarray:
        """
        Get matrix in a specific format.
        
        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        
        Returns
        -------
        out : numpy.ndarray
            Matrix in the desired output format.
        """
        raise NotImplementedError("method is abstract")

    ############## Coancestry/kinship Methods ##############
    def coancestry(self, *args: tuple, **kwargs: dict):
        """
        Retrieve the coancestry between individuals.

        Parameters
        ----------
        args : tuple
            A tuple of matrix indices to access the coancestry.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def kinship(self, *args: tuple, **kwargs: dict):
        """
        Retrieve the kinship between individuals.

        Parameters
        ----------
        args : tuple
            A tuple of matrix indices to access the kinship.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")
    
    def is_positive_semidefinite(
            self, 
            eigvaltol: float
        ) ->  bool:
        """
        Determine whether the coancestry matrix is positive semidefinite.
        
        Parameters
        ----------
        eigvaltol : float
            Eigenvalue tolerance for determining positive semidefiniteness.

        Returns
        -------
        out : bool
            Whether the coancestry matrix is positive semidefinite.
        """
        raise NotImplementedError("method is abstract")

    def apply_jitter(
            self, 
            eigvaltol: float, 
            minjitter: float, 
            maxjitter: float, 
            nattempt: int
        ) -> bool:
        """
        Add a random jitter value to the diagonal of the coancestry matrix until 
        all eigenvalues exceed the provided eigenvalue tolerance.
        This ensures that a matrix can be decomposed using the Cholesky decomposition.
        This routine attempts to apply a jitter 100 times before giving up.

        Parameters
        ----------
        eigvaltol : float
            Eigenvalue tolerance.
        minjitter : float
            Minimum jitter value applied to a diagonal element.
        maxjitter : float
            Maximum jitter value applied to a diagonal element.
        nattempt : int
            Number of jitter application attempts.
        
        Returns
        -------
        out : bool
            Whether the jitter was successfully applied.
        """
        raise NotImplementedError("method is abstract")


    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    @classmethod
    def from_gmat(
            cls, 
            gmat: GenotypeMatrix, 
            **kwargs: dict
        ) -> 'CoancestryMatrix':
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
def is_CoancestryMatrix(v: Any) -> bool:
    """
    Determine whether an object is a CoancestryMatrix.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a CoancestryMatrix object instance.
    """
    return isinstance(v, CoancestryMatrix)

def check_is_CoancestryMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type CoancestryMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, CoancestryMatrix):
        raise TypeError("variable '{0}' must be a CoancestryMatrix".format(vname))
