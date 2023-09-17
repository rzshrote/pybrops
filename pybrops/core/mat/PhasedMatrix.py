"""
Module defining phased matrix interfaces and associated error checking routines.
"""

__all__ = [
    "PhasedMatrix",
    "check_is_PhasedMatrix",
]

from abc import ABCMeta, abstractmethod
from typing import Sequence, Union
import numpy
from numpy.typing import ArrayLike
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.MutableMatrix import MutableMatrix

class PhasedMatrix(MutableMatrix,metaclass=ABCMeta):
    """
    An abstract class for phased matrix wrapper objects.

    A phased matrix is defined as a matrix with a third dimension. This
    interface mostly pertains to phased genotype matrices.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix phase manipulation routines.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############## Phase Metadata Properites ###############
    @property
    @abstractmethod
    def nphase(self) -> int:
        """Number of chromosome phases represented by the matrix."""
        raise NotImplementedError("property is abstract")
    @nphase.setter
    @abstractmethod
    def nphase(self, value: int) -> None:
        """Set number of phases"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def phase_axis(self) -> int:
        """Axis along which phases are stored."""
        raise NotImplementedError("property is abstract")
    @phase_axis.setter
    @abstractmethod
    def phase_axis(self, value: int) -> None:
        """Set phase axis number"""
        raise NotImplementedError("property is abstract")
    
    ############################## Object Methods ##############################

    ######### Matrix element copy-on-manipulation ##########
    @abstractmethod
    def adjoin_phase(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            **kwargs: dict
        ) -> 'PhasedMatrix':
        """
        Add additional elements to the end of the PhasedMatrix along the phase axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhasedMatrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new PhasedMatrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    @abstractmethod
    def delete_phase(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'PhasedMatrix':
        """
        Delete sub-arrays along the phase axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhasedMatrix
            A PhasedMatrix with deleted elements. Note that concat does not occur
            in-place: a new PhasedMatrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    @abstractmethod
    def insert_phase(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            **kwargs: dict
        ) -> 'PhasedMatrix':
        """
        Insert values along the phase axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhasedMatrix
            A PhasedMatrix with values inserted. Note that insert does not occur
            in-place: a new PhasedMatrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    @abstractmethod
    def select_phase(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'PhasedMatrix':
        """
        Select certain values from the Matrix along the phase axis.

        Parameters
        ----------
        indices : ArrayLike (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhasedMatrix
            The output PhasedMatrix with values selected. Note that select does not
            occur in-place: a new PhasedMatrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    @classmethod
    @abstractmethod
    def concat_phase(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'PhasedMatrix':
        """
        Concatenate list of Matrix together along the phase axis.

        Parameters
        ----------
        mats : Sequence of Matrix
            List of PhasedMatrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : PhasedMatrix
            The concatenated PhasedMatrix. Note that concat does not occur in-place:
            a new PhasedMatrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    ######### Matrix element in-place-manipulation #########
    @abstractmethod
    def append_phase(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            **kwargs: dict
        ) -> None:
        """
        Append values to the Matrix along the phase axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def remove_phase(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Remove sub-arrays along the phase axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def incorp_phase(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            **kwargs: dict
        ) -> None:
        """
        Incorporate values along the phase axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_PhasedMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type PhasedMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, PhasedMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,PhasedMatrix.__name__,type(v).__name__))
