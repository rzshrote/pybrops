"""
Module defining interfaces and associated error checking routines for matrices
with trait metadata.
"""

__all__ = [
    "TraitMatrix",
    "check_is_TraitMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from typing import Sequence, Union

import numpy
from numpy.typing import ArrayLike
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.SortableMatrix import SortableMatrix

class TraitMatrix(SortableMatrix,metaclass=ABCMeta):
    """
    An abstract class for matrix wrapper objects with trait metadata.

    The purpose of this abstract class is to provide base functionality for:
        1) trait metadata manipulation routines.
        2) trait manipulation routines.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ###################### Trait data ######################
    @property
    @abstractmethod
    def trait(self) -> object:
        """Trait label."""
        raise NotImplementedError("property is abstract")
    @trait.setter
    @abstractmethod
    def trait(self, value: object) -> None:
        """Set trait label array"""
        raise NotImplementedError("property is abstract")
    
    #################### Trait metadata ####################
    @property
    @abstractmethod
    def ntrait(self) -> int:
        """Number of traits."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def trait_axis(self) -> int:
        """Axis along which traits are stored."""
        raise NotImplementedError("property is abstract")
    
    ############################## Object Methods ##############################

    ######### Matrix element copy-on-manipulation ##########
    @abstractmethod
    def adjoin_trait(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            trait: numpy.ndarray, 
            **kwargs: dict
        ) -> 'TraitMatrix':
        """
        Add additional elements to the end of the Matrix along the trait axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        trait : numpy.ndarray
            Trait names to adjoin to the TraitMatrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : TraitMatrix
            A copy of the TraitMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new TraitMatrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def delete_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'TraitMatrix':
        """
        Delete sub-arrays along the trait axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : TraitMatrix
            A TraitMatrix with deleted elements. Note that concat does not occur
            in-place: a new TraitMatrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def insert_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            trait: numpy.ndarray, 
            **kwargs: dict
        ) -> 'TraitMatrix':
        """
        Insert values along the trait axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        trait : numpy.ndarray
            Trait names to insert into the TraitMatrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : TraitMatrix
            A TraitMatrix with values inserted. Note that insert does not occur
            in-place: a new TraitMatrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def select_trait(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'TraitMatrix':
        """
        Select certain values from the Matrix along the trait axis.

        Parameters
        ----------
        indices : ArrayLike (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : TraitMatrix
            The output TraitMatrix with values selected. Note that select does not
            occur in-place: a new TraitMatrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    @classmethod
    @abstractmethod
    def concat_trait(
            cls,
            mats: Sequence, 
            **kwargs: dict
        ) -> 'TraitMatrix':
        """
        Concatenate list of Matrix together along the trait axis.

        Parameters
        ----------
        mats : Sequence of TraitMatrix
            List of TraitMatrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : TraitMatrix
            The concatenated TraitMatrix. Note that concat does not occur in-place:
            a new TraitMatrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    ######### Matrix element in-place-manipulation #########
    @abstractmethod
    def append_trait(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            trait: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Append values to the TraitMatrix along the trait axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the TraitMatrix.
        trait : numpy.ndarray
            Trait names to append to the TraitMatrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def remove_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Remove sub-arrays along the trait axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def incorp_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            trait: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Incorporate values along the trait axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        trait : numpy.ndarray
            Trait names to incorporate into the Matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ################### Sorting Methods ####################
    @abstractmethod
    def lexsort_trait(
            self, 
            keys: Union[tuple,numpy.ndarray], 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Perform an indirect stable sort using a sequence of keys along the trait
        axis.

        Parameters
        ----------
        keys : A (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        indices : A (N,) ndarray of ints
            Array of indices that sort the keys along the specified axis.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def reorder_trait(
            self, 
            indices: Union[numpy.ndarray,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Reorder elements of the Matrix along the trait axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : A (N,) ndarray of ints, Sequence of ints
            Array of indices that reorder the matrix along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def sort_trait(
            self, 
            keys: Union[tuple,numpy.ndarray], 
            **kwargs: dict
        ) -> None:
        """
        Sort slements of the Matrix along the trait axis using a sequence of
        keys. Note this modifies the Matrix in-place.

        Parameters
        ----------
        keys : A (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_TraitMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type TraitMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, TraitMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,TraitMatrix.__name__,type(v).__name__))
