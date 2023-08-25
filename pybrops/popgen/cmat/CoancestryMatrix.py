"""
Module defining basal coancestry matrix interfaces and associated error checking routines.
"""

__all__ = [
    "CoancestryMatrix",
    "check_is_CoancestryMatrix"
]

from abc import ABCMeta, abstractmethod
from numbers import Real
from typing import Optional, Union
import numpy
from numpy.typing import DTypeLike
from pybrops.core.mat.SquareTaxaMatrix import SquareTaxaMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class CoancestryMatrix(SquareTaxaMatrix,metaclass=ABCMeta):
    """
    An abstract class for coancestry matrices. Coancestry matrices are square.
    Coancestry matrices are related to kinship matrices in the following manner:

    ..math:
        \\mathbf{K} = \\frac{1}{2}\\mathbf{A}

    The purpose of this abstract class is to define base functionality for:
        1) Coancestry matrix value calculation.
        2) Coancestry matrix value access.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############### Genotype Data Properites ###############
    # gmat should be implemented in a sparse version of this matrix.

    ############## Coancestry Data Properites ##############
    # access using mat (inherited from Matrix)

    ############################## Object Methods ##############################

    ################## Matrix conversion ###################
    @abstractmethod
    def mat_asformat(
            self, 
            format: str
        ) -> numpy.ndarray:
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
    @abstractmethod
    def coancestry(
            self, 
            *args: tuple, 
            **kwargs: dict
        ) -> Real:
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

    @abstractmethod
    def kinship(
            self, 
            *args: tuple, 
            **kwargs: dict
        ) -> Real:
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
    
    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
    def max_inbreeding(
            self,
            format: str
        ) -> Real:
        """
        Calculate the maximum attainable inbreeding after one generation for 
        the coancestry matrix. For coancestry, this is equivalent to:
        
        ..math:
            \\max(\\mathrm{trace}(\\mathbf{G}))

        or for kinship, the equivalent is:

        ..math:
            \\max(\\mathrm{trace}(\\mathbf{K}))

        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        
        Returns
        -------
        out : Real
            The maximum attainable inbreeding after one generation.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def min_inbreeding(
            self,
            format: str
        ) -> Real:
        """
        Calculate the minimum attainable inbreeding after one generation for 
        the coancestry matrix. For coancestry, this is equivalent to:
        
        ..math:
            \\frac{1}{\\mathbf{1'G1}}

        or for kinship, the equivalent is:

        ..math:
            \\frac{1}{\\mathbf{1'K1}}

        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        
        Returns
        -------
        out : Real
            The minimum attainable inbreeding after one generation.
        """
        raise NotImplementedError("method is abstract")

    ############## Matrix summary statistics ###############
    @abstractmethod
    def inverse(
            self,
            format: str
        ) -> numpy.ndarray:
        """
        Calculate the inverse of the coancestry matrix.

        Parameters
        ----------
        format : str
            Desired matrix type on which to calculate the inverse. 
            Options are "coancestry", "kinship".

        Returns
        -------
        out : numpy.ndarray
            Inverse of the coancestry or kinship matrix.
        """
        raise NotImplementedError("method is abstract")
    
    @abstractmethod
    def max(
            self,
            format: str,
            axis: Union[int,tuple,None]
    ) -> Union[Real,numpy.ndarray]:
        """
        Calculate the maximum coancestry or kinship for the CoancestryMatrix
        along a specified axis.

        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        axis : int, tuple of ints, None
            Axis along which to find the maximum value.
        
        Returns
        -------
        out : Real, numpy.ndarray
            Maximum coancestry or kinship for the CoancestryMatrix along the 
            specified axis.
        """
        raise NotImplementedError("method is abstract")
    
    @abstractmethod
    def mean(
            self,
            format: str,
            axis: Union[int,tuple,None],
            dtype: Optional[DTypeLike]
        ) -> Real:
        """
        Calculate the mean coancestry or kinship for the CoancestryMatrix
        along a specified axis.

        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        axis : int, tuple of ints, None
            Axis along which to find the mean value.
        dtype : DTypeLike, None
            Type to use in computing the mean. If ``None`` use the native 
            float type.

        Returns
        -------
        out : Real
            Mean coancestry or kinship for the CoancestryMatrix along the 
            specified axis.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def min(
            self,
            format: str,
            axis: Union[int,tuple,None]
    ) -> Union[Real,numpy.ndarray]:
        """
        Calculate the minimum coancestry or kinship for the CoancestryMatrix
        along a specified axis.

        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        axis : int, tuple of ints, None
            Axis along which to find the minimum value.
        
        Returns
        -------
        out : Real, numpy.ndarray
            Minimum coancestry or kinship for the CoancestryMatrix along the 
            specified axis.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################
    @classmethod
    @abstractmethod
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



################################## Utilities ###################################
def check_is_CoancestryMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type CoancestryMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, CoancestryMatrix):
        raise TypeError("variable '{0}' must be a CoancestryMatrix".format(vname))
