"""
Module implementing a dense, scaled matrix and associated error checking routines.
"""

__all__ = [
    "DenseScaledMatrix",
    "check_is_DenseScaledMatrix",
]

import copy
from numbers import Real
from typing import Optional, Union
import numpy
from pybrops.core.error.error_type_numpy import check_is_Real_or_ndarray, check_is_ndarray
from pybrops.core.error.error_type_python import check_is_bool
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len, check_ndarray_ndim
from pybrops.core.mat.DenseMatrix import DenseMatrix
from pybrops.core.mat.ScaledMatrix import ScaledMatrix

class DenseScaledMatrix(
        DenseMatrix,
        ScaledMatrix,
    ):
    """
    A concrete class for dense, scaled matrix wrapper objects.

    The purpose of this concrete class is to provide base functionality for:
        1) Matrix in-place matrix scaling/unscaling routines.
        2) Matrix scaling/unscaling routines
    """

    ########################## Special Object Methods ##########################

    ##################### Constructor ######################
    def __init__(
            self, 
            mat: numpy.ndarray,
            location: Union[numpy.ndarray,Real] = 0.0,
            scale: Union[numpy.ndarray,Real] = 1.0,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseScaledMatrix.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # call DenseMatrix constructor
        super(DenseScaledMatrix, self).__init__(
            mat = mat,
            **kwargs
        )
        # set location and scale
        self.location = location
        self.scale = scale

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseScaledMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseScaledMatrix
            A copy of the DenseScaledMatrix.
        """
        return self.__class__(
            mat = copy.copy(self.mat),
            location = copy.copy(self.location),
            scale = copy.copy(self.scale),
        )
    
    def __deepcopy__(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseScaledMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict, None, default = None
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseScaledMatrix
            A deep copy of the DenseScaledMatrix.
        """
        return self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            location = copy.deepcopy(self.location, memo),
            scale = copy.deepcopy(self.scale, memo),
        )

    ############################ Object Properties #############################

    ############## Scale metadata properties ###############

    @property
    def location(self) -> numpy.ndarray:
        """Location of the matrix data."""
        return self._location
    @location.setter
    def location(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set the location of the matrix data."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "location", 1)
            check_ndarray_axis_len(value, "location", 0, self.mat_shape[-1])
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.mat_shape[-1])
        else:
            check_is_Real_or_ndarray(value, "location")
        self._location = value
    
    @property
    def scale(self) -> numpy.ndarray:
        """Scale of the matrix data."""
        return self._scale
    @scale.setter
    def scale(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set the scale of the matrix data."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "scale", 1)
            check_ndarray_axis_len(value, "scale", 0, self.mat_shape[-1])
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.mat_shape[-1])
        else:
            check_is_Real_or_ndarray(value, "scale")
        self._scale = value

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseScaledMatrix':
        """
        Make a shallow copy of the DenseScaledMatrix.

        Returns
        -------
        out : DenseScaledMatrix
            A shallow copy of the original DenseScaledMatrix.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseScaledMatrix':
        """
        Make a deep copy of the DenseScaledMatrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseScaledMatrix
            A deep copy of the original DenseScaledMatrix.
        """
        return copy.deepcopy(self, memo)

    ################### Scaling methods ####################

    def transform(
            self, 
            mat: numpy.ndarray, 
            copy: bool = False,
        ) -> numpy.ndarray:
        """
        Transform a provided matrix using location and scale parameters within 
        the DenseScaledMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            An array to be centered and scaled.
        copy : bool, default = False
            Whether to copy the input ``mat`` or not.

        Returns
        -------
        out : numpy.ndarray
            A transformed array.
        """
        # type checks
        check_is_ndarray(mat, "mat")
        check_is_bool(copy, "copy")

        # copy matrix if needed
        out = mat.copy() if copy else mat

        # mean center values
        # (?,...,?,n) - (n,) -> (?,...,?,n)
        out -= self.location

        # scale values
        # (?,...,?,n) * (n,) -> (?,...,?,n)
        out *= (1.0 / self.scale)

        return out
    
    def untransform(
            self, 
            mat: numpy.ndarray, 
            copy: bool = False,
        ) -> numpy.ndarray:
        """
        Untransform a provided matrix using location and scale parameters
        within the DenseScaledMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            An array to be unscaled and decentered.
        copy : bool, default = False
            Whether to copy the input ``mat`` or not.

        Returns
        -------
        out : numpy.ndarray
            A transformed array.
        """
        # type checks
        check_is_ndarray(mat, "mat")
        check_is_bool(copy, "copy")

        # copy matrix if needed
        out = mat.copy() if copy else mat

        # unscale values
        # (?,...,?,n) * (n,) -> (?,...,?,n)
        out *= self.scale

        # uncenter values
        # (?,...,?,n) + (n,) -> (?,...,?,n)
        out += self.location

        return out
    
    def rescale(
            self, 
            inplace: bool = True,
        ) -> numpy.ndarray:
        """
        Transform matrix values within the ScaledMatrix to have centered and 
        scaled values.

        Parameters
        ----------
        inplace : bool, default = True
            Whether to modify matrix values in-place.
            If ``True``, then values in DenseScaledMatrix values are scaled internally.
        
        Returns
        -------
        out : numpy.ndarray
            A matrix of centered and scaled values. If ``inplace == True``, then
            a pointer to the raw matrix is returned.
        """
        # type checks
        check_is_bool(inplace, "inplace")

        # get matrix on which to work
        out = self.mat if inplace else self.mat.copy()

        # unscale values
        out *= self.scale
        out += self.location

        # get the axes along which to calculate the location and scale
        axes = tuple(range(out.ndim-1))

        # calculate new location parameters
        # (?,...,?,t) -> (t,)
        new_location = numpy.nanmean(out, axis = axes)

        # calculate new scale parameters
        # (?,...,?,t) -> (t,)
        new_scale = numpy.nanstd(out, axis = axes)
        new_scale[new_scale == 0.0] = 1.0   # if scale == 0.0, set to 1.0 (do not scale)

        # rescale values
        out -= new_location
        out *= (1.0 / new_scale)

        # if in-place modification, update location, scale
        if inplace:
            self.location = new_location    # set location
            self.scale = new_scale          # set scale

        return out
    
    def unscale(
            self, 
            inplace: bool = True,
        ) -> numpy.ndarray:
        """
        Transform matrix values within the ScaledMatrix back to their unscaled 
        and decentered values.

        Parameters
        ----------
        inplace : bool, default = True
            Whether to modify matrix values in-place.
        
        Returns
        -------
        out : numpy.ndarray
            A matrix of unscaled and decentered values. If ``inplace == True``,
            then a pointer to the raw matrix is returned.
        """
        # type checks
        check_is_bool(inplace, "inplace")

        # get matrix on which to work
        out = self.mat if inplace else self.mat.copy()

        # unscale and uncenter values
        # (?,...,?,n) * (n,) -> (?,...,?,n)
        # (?,...,?,n) + (n,) -> (?,...,?,n)
        out *= self.scale
        out += self.location

        if inplace:
            self.scale[:] = 1.0
            self.location[:] = 0.0

        return out

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_DenseScaledMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``DenseScaledMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseScaledMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                DenseScaledMatrix.__name__,
                type(v).__name__
            )
        )
