"""
Module defining interfaces and associated error checking routines for matrices
with variant metadata.
"""

__all__ = [
    "VariantMatrix",
    "check_is_VariantMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from typing import Sequence
from typing import Union

import numpy
from numpy.typing import ArrayLike
from pybrops.core.mat.GroupableMatrix import GroupableMatrix
from pybrops.core.mat.Matrix import Matrix

class VariantMatrix(
        GroupableMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for matrix wrapper objects with variant metadata.

    The purpose of this abstract class is to provide base functionality for:
        1) variant metadata manipulation routines.
        2) variant manipulation routines.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############### Variant Data Properites ################
    @property
    @abstractmethod
    def vrnt_chrgrp(self) -> object:
        """Variant chromosome group label."""
        raise NotImplementedError("property is abstract")
    @vrnt_chrgrp.setter
    @abstractmethod
    def vrnt_chrgrp(self, value: object) -> None:
        """Set variant chromosome group lable array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_phypos(self) -> object:
        """Variant physical position."""
        raise NotImplementedError("property is abstract")
    @vrnt_phypos.setter
    @abstractmethod
    def vrnt_phypos(self, value: object) -> None:
        """Set variant physical position array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_name(self) -> object:
        """Variant name."""
        raise NotImplementedError("property is abstract")
    @vrnt_name.setter
    @abstractmethod
    def vrnt_name(self, value: object) -> None:
        """Set variant name array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_genpos(self) -> object:
        """Variant genetic position."""
        raise NotImplementedError("property is abstract")
    @vrnt_genpos.setter
    @abstractmethod
    def vrnt_genpos(self, value: object) -> None:
        """Set variant genetic position array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_xoprob(self) -> object:
        """Variant crossover sequential probability."""
        raise NotImplementedError("property is abstract")
    @vrnt_xoprob.setter
    @abstractmethod
    def vrnt_xoprob(self, value: object) -> None:
        """Set variant crossover sequential probability array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_hapgrp(self) -> object:
        """Variant haplotype group label."""
        raise NotImplementedError("property is abstract")
    @vrnt_hapgrp.setter
    @abstractmethod
    def vrnt_hapgrp(self, value: object) -> None:
        """Set variant haplotype group label array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_hapalt(self) -> object:
        """Variant haplotype sequence."""
        raise NotImplementedError("property is abstract")
    @vrnt_hapalt.setter
    @abstractmethod
    def vrnt_hapalt(self, value: object) -> None:
        """Set variant haplotype sequence"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_hapref(self) -> object:
        """Variant reference haplotype sequence."""
        raise NotImplementedError("property is abstract")
    @vrnt_hapref.setter
    @abstractmethod
    def vrnt_hapref(self, value: object) -> None:
        """Set variant reference haplotype sequence"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_mask(self) -> object:
        """Variant mask."""
        raise NotImplementedError("property is abstract")
    @vrnt_mask.setter
    @abstractmethod
    def vrnt_mask(self, value: object) -> None:
        """Set variant mask"""
        raise NotImplementedError("property is abstract")
    
    ############# Variant Metadata Properites ##############
    @property
    @abstractmethod
    def nvrnt(self) -> int:
        """Number of variants."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_axis(self) -> int:
        """Axis along which variants are stored."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_chrgrp_name(self) -> object:
        """Variant chromosome group names."""
        raise NotImplementedError("property is abstract")
    @vrnt_chrgrp_name.setter
    @abstractmethod
    def vrnt_chrgrp_name(self, value: object) -> None:
        """Set variant chromosome group name array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_chrgrp_stix(self) -> object:
        """Variant chromosome group start indices."""
        raise NotImplementedError("property is abstract")
    @vrnt_chrgrp_stix.setter
    @abstractmethod
    def vrnt_chrgrp_stix(self, value: object) -> None:
        """Set variant chromosome group start indices array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_chrgrp_spix(self) -> object:
        """Variant chromosome group stop indices."""
        raise NotImplementedError("property is abstract")
    @vrnt_chrgrp_spix.setter
    @abstractmethod
    def vrnt_chrgrp_spix(self, value: object) -> None:
        """Set variant chromosome group stop indices array"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def vrnt_chrgrp_len(self) -> object:
        """Variant chromosome group length."""
        raise NotImplementedError("property is abstract")
    @vrnt_chrgrp_len.setter
    @abstractmethod
    def vrnt_chrgrp_len(self, value: object) -> None:
        """Set variant chromosome group length array"""
        raise NotImplementedError("property is abstract")
    
    ############################## Object Methods ##############################

    ######### Matrix element copy-on-manipulation ##########
    @abstractmethod
    def adjoin_vrnt(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            vrnt_name: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            vrnt_xoprob: numpy.ndarray, 
            vrnt_hapgrp: numpy.ndarray, 
            vrnt_mask: numpy.ndarray, 
            **kwargs: dict
        ) -> 'VariantMatrix':
        """
        Add additional elements to the end of the VariantMatrix along the variant axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to adjoin to the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to adjoin to the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to adjoin to the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to adjoin to the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to adjoin to the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to adjoin to the Matrix.
        vrnt_mask : numpy.ndarray
            Variant mask to adjoin to the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : VariantMatrix
            A copy of the VariantMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new VariantMatrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    @abstractmethod
    def delete_vrnt(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'VariantMatrix':
        """
        Delete sub-arrays along the variant axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : VariantMatrix
            A VariantMatrix with deleted elements. Note that delete does not occur
            in-place: a new VariantMatrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    @abstractmethod
    def insert_vrnt(
            self, 
            obj: Union[int,slice,Sequence], 
            values: ArrayLike, 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            vrnt_name: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            vrnt_xoprob: numpy.ndarray, 
            vrnt_hapgrp: numpy.ndarray, 
            vrnt_mask: numpy.ndarray, 
            **kwargs: dict
        ) -> 'VariantMatrix':
        """
        Insert values along the variant axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : array_like
            Values to insert into the matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to insert into the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to insert into the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to insert into the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to insert into the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to insert into the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to insert into the Matrix.
        vrnt_mask : numpy.ndarray
            Variant mask to insert into the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : VariantMatrix
            A VariantMatrix with values inserted. Note that insert does not occur
            in-place: a new VariantMatrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    @abstractmethod
    def select_vrnt(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'VariantMatrix':
        """
        Select certain values from the VariantMatrix along the variant axis.

        Parameters
        ----------
        indices : ArrayLike (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : VariantMatrix
            The output VariantMatrix with values selected. Note that select does not
            occur in-place: a new VariantMatrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    @classmethod
    @abstractmethod
    def concat_vrnt(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'VariantMatrix':
        """
        Concatenate list of VariantMatrix together along the variant axis.

        Parameters
        ----------
        mats : Sequence of VariantMatrix
            List of VariantMatrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : VariantMatrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new VariantMatrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    ######### Matrix element in-place-manipulation #########
    @abstractmethod
    def append_vrnt(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            vrnt_name: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            vrnt_xoprob: numpy.ndarray, 
            vrnt_hapgrp: numpy.ndarray, 
            vrnt_mask: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Append values to the VariantMatrix along the variant axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the VariantMatrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to append to the VariantMatrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to append to the VariantMatrix.
        vrnt_name : numpy.ndarray
            Variant names to append to the VariantMatrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to append to the VariantMatrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to append to the VariantMatrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to append to the VariantMatrix.
        vrnt_mask : numpy.ndarray
            Variant mask to append to the VariantMatrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def remove_vrnt(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Remove sub-arrays along the variant axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def incorp_vrnt(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            vrnt_chrgrp: numpy.ndarray, 
            vrnt_phypos: numpy.ndarray, 
            vrnt_name: numpy.ndarray, 
            vrnt_genpos: numpy.ndarray, 
            vrnt_xoprob: numpy.ndarray, 
            vrnt_hapgrp: numpy.ndarray, 
            vrnt_mask: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Incorporate values along the variant axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the VariantMatrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to incorporate into the VariantMatrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to incorporate into the VariantMatrix.
        vrnt_name : numpy.ndarray
            Variant names to incorporate into the VariantMatrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to incorporate into the VariantMatrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to incorporate into the VariantMatrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to incorporate into the VariantMatrix.
        vrnt_mask : numpy.ndarray
            Variant mask to incorporate into the VariantMatrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ################### Sorting Methods ####################
    @abstractmethod
    def lexsort_vrnt(
            self, 
            keys: Union[tuple,numpy.ndarray], 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Perform an indirect stable sort using a sequence of keys along the
        variant axis.

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
    def reorder_vrnt(
            self, 
            indices: Union[numpy.ndarray,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Reorder elements of the Matrix along the variant axis using an array of
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
    def sort_vrnt(
            self, 
            keys: Union[tuple,numpy.ndarray], 
            **kwargs: dict
        ) -> None:
        """
        Sort slements of the Matrix along the variant axis using a sequence of
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
    def group_vrnt(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Sort the VariantMatrix along the variant axis, then populate grouping 
        indices for the variant axis.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def ungroup_vrnt(
            self,
            **kwargs: dict
        ) -> None:
        """
        Ungroup the VariantMatrix along the variant axis by removing variant 
        group metadata.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def is_grouped_vrnt(
            self, 
            **kwargs: dict
        ) -> bool:
        """
        Determine whether the Matrix has been sorted and grouped along the
        variant axis.

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
def check_is_VariantMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type VariantMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, VariantMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,VariantMatrix.__name__,type(v).__name__))
