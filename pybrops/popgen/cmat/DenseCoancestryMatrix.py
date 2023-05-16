"""
Module providing dense coancestry matrix implementations and associated error checking routines.
"""

from numbers import Real
from typing import Any, Optional, Union
import numpy
import warnings
from numpy.typing import DTypeLike
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_python import check_all_equal
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_numpy import check_ndarray_dtype
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_object
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_value_python import check_str_value
from pybrops.core.mat.DenseSquareTaxaMatrix import DenseSquareTaxaMatrix
from pybrops.popgen.cmat.CoancestryMatrix import CoancestryMatrix

class DenseCoancestryMatrix(DenseSquareTaxaMatrix,CoancestryMatrix):
    """
    A concrete class for dense coancestry matrices. Coancestry matrices are square.

    The purpose of this concrete class is to implement functionality for:
        1) Dense coancestry matrix value calculation.
        2) Dense coancestry matrix value access.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseCoancestryMatrix class.

        Parameters
        ----------
        mat : numpy.ndarray
            Coancestry matrix of shape ``(n,n)``.

            Where:

            - ``n`` is the number of taxa.
        taxa : numpy.ndarray, None
        taxa_grp : numpy.ndarray, None
        kwargs : dict
            Additional keyword arguments.
        """
        # call constructor for DenseTaxaMatrix
        super(DenseCoancestryMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Coancestry Data Properites ##############
    @DenseSquareTaxaMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        check_is_ndarray(value, "mat")
        check_all_equal(value.shape, "mat.shape")
        check_ndarray_dtype(value, "mat", numpy.float64)
        check_ndarray_ndim(value, "mat", 2)
        self._mat = value

    ################# Taxa Data Properites #################
    @DenseSquareTaxaMatrix.taxa.setter
    def taxa(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa label array"""
        if value is not None:
            check_is_ndarray(value, "taxa")
            check_ndarray_dtype_is_object(value, "taxa")
            check_ndarray_ndim(value, "taxa", 1)
            check_ndarray_axis_len(value, "taxa", 0, self._mat.shape[0])
        self._taxa = value

    @DenseSquareTaxaMatrix.taxa_grp.setter
    def taxa_grp(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group label array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp")
            check_ndarray_dtype(value, "taxa_grp", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp", 1)
            check_ndarray_axis_len(value, "taxa_grp", 0, self._mat.shape[0])
        self._taxa_grp = value

    ############### Taxa Metadata Properites ###############
    @DenseSquareTaxaMatrix.taxa_grp_name.setter
    def taxa_grp_name(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group name array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp_name")
            check_ndarray_dtype(value, "taxa_grp_name", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_name", 1)
        self._taxa_grp_name = value

    @DenseSquareTaxaMatrix.taxa_grp_stix.setter
    def taxa_grp_stix(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group start indices array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp_stix")
            check_ndarray_dtype(value, "taxa_grp_stix", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_stix", 1)
        self._taxa_grp_stix = value

    @DenseSquareTaxaMatrix.taxa_grp_spix.setter
    def taxa_grp_spix(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group stop indices array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp_spix")
            check_ndarray_dtype(value, "taxa_grp_spix", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_spix", 1)
        self._taxa_grp_spix = value

    @DenseSquareTaxaMatrix.taxa_grp_len.setter
    def taxa_grp_len(self, value: Union[numpy.ndarray,None]) -> None:
        if value is not None:
            check_is_ndarray(value, "taxa_grp_len")
            check_ndarray_dtype(value, "taxa_grp_len", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_len", 1)
        self._taxa_grp_len = value

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
        # check values
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")

        # process string
        if format == "coancestry":
            return self._mat.copy()
        elif format == "kinship":
            return 0.5 * self._mat

    ############## Coancestry/kinship Methods ##############
    def coancestry(self, *args, **kwargs: dict):
        """
        Retrieve the coancestry between individuals.

        Parameters
        ----------
        args : tuple
            A tuple of matrix indices to access the coancestry.
        kwargs : dict
            Additional keyword arguments.
        """
        return self._mat[args]

    def kinship(self, *args, **kwargs: dict):
        """
        Retrieve the kinship between individuals.

        Parameters
        ----------
        args : tuple
            A tuple of matrix indices to access the kinship.
        kwargs : dict
            Additional keyword arguments.
        """
        return 0.5 * self._mat[args]

    def is_positive_semidefinite(
            self, 
            eigvaltol = 2e-14
        ) -> bool:
        """
        Determine whether the coancestry matrix is positive semidefinite.
        
        Parameters
        ----------
        eigvaltol : float
            Eigenvalue tolerance for determining positive semidefiniteness.
            If provided eigenvalue tolerance is less than zero, the tolerance
            is set to 0.0.

        Returns
        -------
        out : bool
            Whether the coancestry matrix is positive semidefinite.
        """
        # if tolerance is less than zero, set to zero
        if eigvaltol < 0.0:
            eigvaltol = 0.0
        return numpy.all(numpy.linalg.eigvals(self._mat) >= eigvaltol)
    
    def apply_jitter(
            self, 
            eigvaltol = 2e-14, 
            minjitter = 1e-10, 
            maxjitter = 1e-6, 
            nattempt = 100
        ) -> bool:
        """
        Add a random jitter value to the diagonal of the coancestry matrix until 
        all eigenvalues exceed the provided eigenvalue tolerance.
        This ensures that a matrix can be decomposed using the Cholesky decomposition.
        This routine attempts to apply a jitter 100 times before giving up.

        Parameters
        ----------
        eigvaltol : float
            Eigenvalue tolerance for determining positive semidefiniteness.
            If provided eigenvalue tolerance is less than zero, the tolerance
            is set to 0.0.
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
        diagix = numpy.diag_indices_from(self._mat)
        mat_diag_old = self._mat[diagix].copy()
        counter = 0
        is_not_posdef = not self.is_positive_semidefinite(eigvaltol)

        # attempt to apply jitter in nattempt or less attempts
        while (is_not_posdef) and (counter < nattempt):
            self._mat[diagix] = mat_diag_old + numpy.random.uniform(minjitter, maxjitter, len(mat_diag_old))
            is_not_posdef = not self.is_positive_semidefinite(eigvaltol)
            counter += 1
        
        # if we still weren't able to find appropriate jitter, then give old diagonals and warn
        if is_not_posdef:
            self._mat[diagix] = mat_diag_old
            warnings.warn("Unable to successfully apply jitter to meet eigenvalue tolerance")
            return False
        
        return True

    def max_inbreeding(
            self,
            format: str
        ) -> Real:
        """
        Calculate the maximum attainable inbreeding after one generation for 
        the coancestry matrix. For coancestry, this is equivalent to:
        
        ..math:
            \\max(\\mathrm{diag}(\\mathbf{G}))

        or for kinship, the equivalent is:

        ..math:
            \\max(\\mathrm{diag}(\\mathbf{K}))

        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        
        Returns
        -------
        out : Real
            The maximum attainable inbreeding after one generation.
        """
        # check values
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")

        # get output in coancestry format
        out = self.mat.diagonal().max()

        # if requested format is kinship, then convert coancestry to kinship
        if format == "kinship":
            out = 0.5 * out
        
        return out

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
        # check values
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")

        # calculate G inverse
        Ginv = numpy.linalg.inv(self.mat)

        # calculate min inbreeding
        out = 1.0 / Ginv.sum()

        # if requested format is kinship, then convert coancestry to kinship
        if format == "kinship":
            out = 0.5 * out
        
        return out

    ############## Matrix summary statistics ###############
    def inverse(
            self,
            format: str = "coancestry"
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
        # type checks
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")

        # invert matrix
        if format == "coancestry":
            out = numpy.linalg.inv(self._mat)
        elif format == "kinship":
            out = numpy.linalg.inv(0.5 * self._mat)
        
        return out

    def max(
            self,
            format: str = "coancestry",
            axis: Union[int,tuple,None] = None
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
        # type checks
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")
        
        # calculate mean values
        out = self._mat.max(axis = axis)

        # process string and apply transformations
        if format == "kinship":
            out *= 0.5
        
        return out

    def mean(
            self,
            format: str = "coancestry",
            axis: Union[int,tuple,None] = None,
            dtype: Optional[DTypeLike] = None
        ) -> Real:
        """
        Calculate the mean coancestry or kinship for the CoancestryMatrix.

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
            Mean coancestry or kinship for the CoancestryMatrix.
        """
        # type checks
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")
        
        # calculate mean values
        out = self._mat.mean(axis = axis, dtype = dtype)

        # process string and apply transformations
        if format == "kinship":
            out *= 0.5
        
        return out

    def min(
            self,
            format: str = "coancestry",
            axis: Union[int,tuple,None] = None
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
        # type checks
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")
        
        # calculate mean values
        out = self._mat.min(axis = axis)

        # process string and apply transformations
        if format == "kinship":
            out *= 0.5
        
        return out



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseCoancestryMatrix(v: Any) -> bool:
    """
    Determine whether an object is a DenseCoancestryMatrix.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseCoancestryMatrix object instance.
    """
    return isinstance(v, DenseCoancestryMatrix)

def check_is_DenseCoancestryMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type DenseCoancestryMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseCoancestryMatrix):
        raise TypeError("variable '{0}' must be a DenseCoancestryMatrix".format(vname))
