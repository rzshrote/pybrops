"""
Module defining interfaces and associated error checking routines for matrices
with taxa metadata.
"""

__all__ = [
    "TaxaMatrix",
    "check_is_TaxaMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from typing import Sequence
from typing import Union

import numpy
from numpy.typing import ArrayLike
from pybrops.core.mat.GroupableMatrix import GroupableMatrix
from pybrops.core.mat.Matrix import Matrix

class TaxaMatrix(GroupableMatrix,metaclass=ABCMeta):
    """
    An abstract class for matrix wrapper objects with taxa metadata.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix taxa metadata.
        2) Matrix taxa routines.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ################# Taxa Data Properites #################
    @property
    @abstractmethod
    def taxa(self) -> object:
        """Taxa label."""
        raise NotImplementedError("property is abstract")
    @taxa.setter
    @abstractmethod
    def taxa(self, value: object) -> None:
        """Set taxa label array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def taxa_grp(self) -> object:
        """Taxa group label."""
        raise NotImplementedError("property is abstract")
    @taxa_grp.setter
    @abstractmethod
    def taxa_grp(self, value: object) -> None:
        """Set taxa group label array"""
        raise NotImplementedError("property is abstract")
    
    ############### Taxa Metadata Properites ###############
    @property
    @abstractmethod
    def ntaxa(self) -> int:
        """Number of taxa."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def taxa_axis(self) -> int:
        """Axis along which taxa are stored."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def taxa_grp_name(self) -> object:
        """Taxa group name."""
        raise NotImplementedError("property is abstract")
    @taxa_grp_name.setter
    @abstractmethod
    def taxa_grp_name(self, value: object) -> None:
        """Set taxa group name array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def taxa_grp_stix(self) -> object:
        """Taxa group start index."""
        raise NotImplementedError("property is abstract")
    @taxa_grp_stix.setter
    @abstractmethod
    def taxa_grp_stix(self, value: object) -> None:
        """Set taxa group start indices array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def taxa_grp_spix(self) -> object:
        """Taxa group stop index."""
        raise NotImplementedError("property is abstract")
    @taxa_grp_spix.setter
    @abstractmethod
    def taxa_grp_spix(self, value: object) -> None:
        """Set taxa group stop indices array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def taxa_grp_len(self) -> object:
        """Taxa group length."""
        raise NotImplementedError("property is abstract")
    @taxa_grp_len.setter
    @abstractmethod
    def taxa_grp_len(self, value: object) -> None:
        """Set taxa group length array"""
        raise NotImplementedError("property is abstract")
    
    ############################## Object Methods ##############################

    ######### Matrix element copy-on-manipulation ##########
    @abstractmethod
    def adjoin_taxa(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            taxa: numpy.ndarray, 
            taxa_grp: numpy.ndarray, 
            **kwargs: dict
        ) -> 'TaxaMatrix':
        """
        Add additional elements to the end of the Matrix along the taxa axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : TaxaMatrix
            A copy of TaxaMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new TaxaMatrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def delete_taxa(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'TaxaMatrix':
        """
        Delete sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : TaxaMatrix
            A TaxaMatrix with deleted elements. Note that concat does not occur
            in-place: a new TaxaMatrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def insert_taxa(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            taxa: numpy.ndarray, 
            taxa_grp: numpy.ndarray, 
            **kwargs: dict
        ) -> 'TaxaMatrix':
        """
        Insert values along the taxa axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        taxa : numpy.ndarray
            Taxa names to insert into the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : TaxaMatrix
            A TaxaMatrix with values inserted. Note that insert does not occur
            in-place: a new TaxaMatrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def select_taxa(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'TaxaMatrix':
        """
        Select certain values from the Matrix along the taxa axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : TaxaMatrix
            The output TaxaMatrix with values selected. Note that select does not
            occur in-place: a new TaxaMatrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    @classmethod
    @abstractmethod
    def concat_taxa(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'TaxaMatrix':
        """
        Concatenate list of Matrix together along the taxa axis.

        Parameters
        ----------
        mats : Sequence of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : TaxaMatrix
            The concatenated TaxaMatrix. Note that concat does not occur in-place:
            a new TaxaMatrix is allocated and filled.
        """
        raise NotImplementedError("class method is abstract")

    ######### Matrix element in-place-manipulation #########
    @abstractmethod
    def append_taxa(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            taxa: numpy.ndarray, 
            taxa_grp: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Append values to the Matrix along the taxa axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        taxa : numpy.ndarray
            Taxa names to append to the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to append to the Matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def remove_taxa(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Remove sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def incorp_taxa(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            taxa: numpy.ndarray, 
            taxa_grp: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Incorporate values along the taxa axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        taxa : numpy.ndarray
            Taxa names to incorporate into the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to incorporate into the Matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ################### Sorting Methods ####################
    @abstractmethod
    def lexsort_taxa(
            self, 
            keys: Union[tuple,numpy.ndarray], 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Perform an indirect stable sort using a sequence of keys along the taxa
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
    def reorder_taxa(
            self, 
            indices: Union[numpy.ndarray,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Reorder elements of the Matrix along the taxa axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : A (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def sort_taxa(
            self, 
            keys: Union[tuple,numpy.ndarray], 
            **kwargs: dict
        ) -> None:
        """
        Sort slements of the Matrix along the taxa axis using a sequence of
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

    ################### Grouping Methods ###################
    @abstractmethod
    def group_taxa(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Sort the Matrix along the taxa axis, then populate grouping indices for
        the taxa axis.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def ungroup_taxa(
            self,
            **kwargs: dict
        ) -> None:
        """
        Ungroup the TaxaMatrix along the taxa axis by removing taxa group 
        metadata.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def is_grouped_taxa(
            self, 
            **kwargs: dict
        ) -> bool:
        """
        Determine whether the Matrix has been sorted and grouped along the taxa
        axis.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        grouped : bool
            True or False indicating whether the Matrix has been sorted and
            grouped.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_TaxaMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type TaxaMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, TaxaMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,TaxaMatrix.__name__,type(v).__name__))
