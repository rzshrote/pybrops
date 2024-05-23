"""
Module defining interfaces and associated error checking routines for matrices
with axes that are scaled.
"""

__all__ = [
    "ScaledMatrix",
    "check_is_ScaledMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from numbers import Real
from typing import Union
import numpy
from pybrops.core.mat.Matrix import Matrix

class ScaledMatrix(
        Matrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for scaled matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix in-place matrix scaling/unscaling routines.
        2) Matrix scaling/unscaling routines
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############## Scale metadata properties ###############
    @property
    @abstractmethod
    def location(self) -> numpy.ndarray:
        """Location of the matrix data."""
        raise NotImplementedError("property is abstract")
    @location.setter
    @abstractmethod
    def location(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set the location of the matrix data."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def scale(self) -> numpy.ndarray:
        """Scale of the matrix data."""
        raise NotImplementedError("property is abstract")
    @scale.setter
    @abstractmethod
    def scale(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set the scale of the matrix data."""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    ################### Scaling methods ####################
    @abstractmethod
    def transform(
            self,
            mat: numpy.ndarray,
            copy: bool,
        ) -> numpy.ndarray:
        """
        Transform a provided matrix using location and scale parameters within 
        the ScaledMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            An array to be centered and scaled.
        copy : bool
            Whether to copy the input ``mat`` or not.

        Returns
        -------
        out : numpy.ndarray
            A transformed array.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def untransform(
            self,
            mat: numpy.ndarray,
            copy: bool,
        ) -> numpy.ndarray:
        """
        Untransform a provided matrix using location and scale parameters
        within the ScaledMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            An array to be unscaled and decentered.
        copy : bool
            Whether to copy the input ``mat`` or not.

        Returns
        -------
        out : numpy.ndarray
            A transformed array.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def rescale(
            self, 
            inplace: bool,
        ) -> numpy.ndarray:
        """
        Transform matrix values within the ScaledMatrix to have centered and 
        scaled values.

        Parameters
        ----------
        inplace : bool
            Whether to modify matrix values in-place.
        
        Returns
        -------
        out : numpy.ndarray
            A matrix of centered and scaled values. If ``inplace == True``, then
            a pointer to the raw matrix is returned.
        """
        raise NotImplementedError("method is abstract")
    
    @abstractmethod
    def unscale(
            self,
            inplace: bool,
        ) -> numpy.ndarray:
        """
        Transform matrix values within the ScaledMatrix back to their unscaled 
        and decentered values.

        Parameters
        ----------
        inplace : bool
            Whether to modify matrix values in-place.
        
        Returns
        -------
        out : numpy.ndarray
            A matrix of unscaled and decentered values. If ``inplace == True``,
            then a pointer to the raw matrix is returned.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_ScaledMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``ScaledMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, ScaledMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                ScaledMatrix.__name__,
                type(v).__name__
            )
        )
