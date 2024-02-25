"""
Module defining basal Matrix interfaces and associated error checking routines.
"""

__all__ = [
    "Matrix",
    "check_is_Matrix",
]

from abc import ABCMeta, abstractmethod
import numpy
from typing import Iterator, Union

from typing import Sequence
from numpy.typing import ArrayLike

from pybrops.core.io.HDF5InputOutput import HDF5InputOutput

class Matrix(
        HDF5InputOutput,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for matrix wrapper objects.

    The purpose of this abstract class is to define base functionality for:
    
        1. Matrix mathematical operators
        2. Matrix logical & bitwise operators
        3. Matrix container operators
        4. Matrix copy operators
        5. Matrix read-only matrix shape changing routines.

    The shape of a Matrix should be immutable.
    """

    ########################## Special Object Methods ##########################

    ############## Forward numeric operators ###############
    @abstractmethod
    def __add__(self, value: object) -> object:
        """
        Elementwise add matrices
        
        Parameters
        ----------
        value : object
            Object which to add.
        
        Returns
        -------
        out : object
            An object resulting from the addition.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __sub__(self, value: object) -> object:
        """
        Elementwise subtract matrices
        
        Parameters
        ----------
        value : object
            Object which to subtract.
        
        Returns
        -------
        out : object
            An object resulting from the subtraction.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __mul__(self, value: object) -> object:
        """
        Elementwise multiply matrices
        
        Parameters
        ----------
        value : object
            Object which to multiply.
        
        Returns
        -------
        out : object
            An object resulting from the multiplication.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __matmul__(self, value: object) -> object:
        """
        Multiply matrices
        
        Parameters
        ----------
        value : object
            Object which to matrix multiply.
        
        Returns
        -------
        out : object
            An object resulting from the matrix multiplication.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __truediv__(self, value: object) -> object:
        """
        Elementwise divide matrices
        
        Parameters
        ----------
        value : object
            Object which to divide.
        
        Returns
        -------
        out : object
            An object resulting from the division.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __floordiv__(self, value: object) -> object:
        """
        Elementwise floor divide matrices
        
        Parameters
        ----------
        value : object
            Object which to floor divide.
        
        Returns
        -------
        out : object
            An object resulting from the floor division.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __mod__(self, value: object) -> object:
        """
        Elementwise modulo matrices

        Parameters
        ----------
        value : object
            Object which to modulo.
        
        Returns
        -------
        out : object
            An object resulting from the modulo.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __divmod__(self, value: object) -> object:
        """
        Elementwise divmod matrices
        
        Parameters
        ----------
        value : object
            Object which to divmod.
        
        Returns
        -------
        out : object
            An object resulting from the divmod.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __pow__(self, value: object) -> object:
        """
        Elementwise exponent matrices
        
        Parameters
        ----------
        value : object
            Object which to exponentiate.
        
        Returns
        -------
        out : object
            An object resulting from the exponentiation.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __lshift__(self, value: object) -> object:
        """
        Elementwise bitwise left shift matrices
        
        Parameters
        ----------
        value : object
            Object which to left shift.
        
        Returns
        -------
        out : object
            An object resulting from the left shift.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rshift__(self, value: object) -> object:
        """
        Elementwise bitwise right shift matrices
        
        Parameters
        ----------
        value : object
            Object which to right shift.
        
        Returns
        -------
        out : object
            An object resulting from the right shift.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __and__(self, value: object) -> object:
        """
        Elementwise bitwise AND matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise AND.
        
        Returns
        -------
        out : object
            An object resulting from the bitwise AND.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __xor__(self, value: object) -> object:
        """
        Elementwise bitwise XOR matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise XOR.
        
        Returns
        -------
        out : object
            An object resulting from the bitwise XOR.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __or__(self, value: object) -> object:
        """
        Elementwise bitwise OR matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise OR.
        
        Returns
        -------
        out : object
            An object resulting from the bitwise OR.
        """
        raise NotImplementedError("method is abstract")

    ############# Backwards numeric operators ##############
    @abstractmethod
    def __radd__(self, value: object) -> object:
        """
        Reverse elementwise add matrices
        
        Parameters
        ----------
        value : object
            Object which to add.
        
        Returns
        -------
        out : object
            An object resulting from the addition.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rsub__(self, value: object) -> object:
        """
        Reverse elementwise subtract matrices
        
        Parameters
        ----------
        value : object
            Object which to subtract.
        
        Returns
        -------
        out : object
            An object resulting from the subtraction.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rmul__(self, value: object) -> object:
        """
        Reverse elementwise multiply matrices
        
        Parameters
        ----------
        value : object
            Object which to multiply.
        
        Returns
        -------
        out : object
            An object resulting from the multiplication.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rmatmul__(self, value: object) -> object:
        """
        Reverse multiply matrices
        
        Parameters
        ----------
        value : object
            Object which to matrix multiply.
        
        Returns
        -------
        out : object
            An object resulting from the matrix multiplication.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rtruediv__(self, value: object) -> object:
        """
        Reverse elementwise divide matrices
        
        Parameters
        ----------
        value : object
            Object which to divide.
        
        Returns
        -------
        out : object
            An object resulting from the division.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rfloordiv__(self, value: object) -> object:
        """
        Reverse elementwise floor divide matrices
        
        Parameters
        ----------
        value : object
            Object which to floor divide.
        
        Returns
        -------
        out : object
            An object resulting from the floor division.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rmod__(self, value: object) -> object:
        """
        Reverse elementwise modulus matrices
        
        Parameters
        ----------
        value : object
            Object which to modulo.
        
        Returns
        -------
        out : object
            An object resulting from the modulo.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rdivmod__(self, value: object) -> object:
        """
        Reverse elementwise divmod matrices
        
        Parameters
        ----------
        value : object
            Object which to divmod.
        
        Returns
        -------
        out : object
            An object resulting from the divmod.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rlshift__(self, value: object) -> object:
        """
        Reverse elementwise bitwise left shift matrices
        
        Parameters
        ----------
        value : object
            Object which to left shift.
        
        Returns
        -------
        out : object
            An object resulting from the left shift.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rrshift__(self, value: object) -> object:
        """
        Reverse elementwise bitwise right shift matrices
        
        Parameters
        ----------
        value : object
            Object which to right shift.
        
        Returns
        -------
        out : object
            An object resulting from the right shift.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rand__(self, value: object) -> object:
        """
        Reverse elementwise bitwise AND matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise AND.
        
        Returns
        -------
        out : object
            An object resulting from the bitwise AND.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __rxor__(self, value: object) -> object:
        """
        Reverse elementwise bitwise XOR matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise XOR.
        
        Returns
        -------
        out : object
            An object resulting from the bitwise XOR.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __ror__(self, value: object) -> object:
        """
        Reverse elementwise bitwise OR matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise OR.
        
        Returns
        -------
        out : object
            An object resulting from the bitwise OR.
        """
        raise NotImplementedError("method is abstract")

    ############# Augmented numeric operators ##############
    @abstractmethod
    def __iadd__(self, value: object) -> None:
        """
        Elementwise add assign matrices
        
        Parameters
        ----------
        value : object
            Object which to add.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __isub__(self, value: object) -> None:
        """
        Elementwise subtract assign matrices
        
        Parameters
        ----------
        value : object
            Object which to subtract.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __imul__(self, value: object) -> None:
        """
        Elementwise multiply assign matrices
        
        Parameters
        ----------
        value : object
            Object which to multiply.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __imatmul__(self, value: object) -> None:
        """
        Multiply assign matrices
        
        Parameters
        ----------
        value : object
            Object which to matrix multiply.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __itruediv__(self, value: object) -> None:
        """
        Elementwise true divide assign matrices
        
        Parameters
        ----------
        value : object
            Object which to divide.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __ifloordiv__(self, value: object) -> None:
        """
        Elementwise floor divide assign matrices
        
        Parameters
        ----------
        value : object
            Object which to floor divide.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __imod__(self, value: object) -> None:
        """
        Elementwise modulus assign matrices
        
        Parameters
        ----------
        value : object
            Object which to modulo.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __ipow__(self, value: object) -> None:
        """
        Elementwise exponent assign matrices
        
        Parameters
        ----------
        value : object
            Object which to exponentiate.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __ilshift__(self, value: object) -> None:
        """
        Elementwise left shift assign matrices
        
        Parameters
        ----------
        value : object
            Object which to left shift.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __irshift__(self, value: object) -> None:
        """
        Elementwise right shift assign matrices
        
        Parameters
        ----------
        value : object
            Object which to right shift.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __iand__(self, value: object) -> None:
        """
        Elementwise bitwise AND assign matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise AND.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __ixor__(self, value: object) -> None:
        """
        Elementwise bitwise XOR assign matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise XOR.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __ior__(self, value: object) -> None:
        """
        Elementwise bitwise OR assign matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise OR.
        """
        raise NotImplementedError("method is abstract")

    ################## Logical operators ###################
    @abstractmethod
    def __lt__(self, value: object) -> object:
        """
        Elementwise less than comparison matrices
        
        Parameters
        ----------
        value : object
            Object which to compare.
        
        Returns
        -------
        out : object
            An object resulting from the comparision.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __le__(self, value: object) -> object:
        """
        Elementwise less than or equal to comparison matrices
        
        Parameters
        ----------
        value : object
            Object which to compare.
        
        Returns
        -------
        out : object
            An object resulting from the comparision.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __eq__(self, value: object) -> object:
        """
        Elementwise equal to comparison matrices
        
        Parameters
        ----------
        value : object
            Object which to compare.
        
        Returns
        -------
        out : object
            An object resulting from the comparision.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __ne__(self, value: object) -> object:
        """
        Elementwise not equal to comparison matrices
        
        Parameters
        ----------
        value : object
            Object which to compare.
        
        Returns
        -------
        out : object
            An object resulting from the comparision.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __gt__(self, value: object) -> object:
        """
        Elementwise greater than comparison matrices
        
        Parameters
        ----------
        value : object
            Object which to compare.
        
        Returns
        -------
        out : object
            An object resulting from the comparision.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __ge__(self, value: object) -> object:
        """
        Elementwise greater than or equal to comparison matrices
        
        Parameters
        ----------
        value : object
            Object which to compare.
        
        Returns
        -------
        out : object
            An object resulting from the comparision.
        """
        raise NotImplementedError("method is abstract")

    ################# Container operators ##################
    @abstractmethod
    def __len__(self) -> int:
        """
        Get the length of the Matrix.

        Returns
        -------
        out : int
            The length of the Matrix.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __getitem__(self, key: object) -> object:
        """
        Get an item from the Matrix.

        Parameters
        ----------
        key : object
            Key of the item which to get.
        
        Returns
        -------
        out : object
            Item at the provided key.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __setitem__(self, key: object, value: object) -> None:
        """
        Set an item in the Matrix.

        Parameters
        ----------
        key : object
            Key of the item which to set.
        value : object
            Value of the item to wich to set.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __delitem__(self, key: object) -> None:
        """
        Delete an item in the Matrix.

        Parameters
        ----------
        key : object
            Key of the item which to delete.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __iter__(self) -> Iterator:
        """
        Get an iterator for the Matrix.

        Returns
        -------
        out : Iterator
            An iterator for the Matrix.
        """
        raise NotImplementedError("method is abstract")

    #################### Matrix copying ####################
    @abstractmethod
    def __copy__(
            self
        ) -> 'Matrix':
        """
        Make a shallow copy of the Matrix.

        Returns
        -------
        out : Matrix
            A shallow copy of the original Matrix.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'Matrix':
        """
        Make a deep copy of the Matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : Matrix
            A deep copy of the original Matrix.
        """
        raise NotImplementedError("method is abstract")

    ########### Miscellaneous special functions ############
    # @abstractmethod
    # def __repr__(self) -> str:
    #     """Return repr(self)."""
    #     raise NotImplementedError("method is abstract")

    ############################ Object Properties #############################
    @property
    @abstractmethod
    def mat(self) -> object:
        """Pointer to raw matrix object."""
        raise NotImplementedError("property is abstract")
    @mat.setter
    @abstractmethod
    def mat(self, value: object) -> None:
        """Set pointer to raw matrix object"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def mat_ndim(self) -> int:
        """Number of dimensions of the raw matrix."""
        raise NotImplementedError("property is abstract")
    @mat_ndim.setter
    @abstractmethod
    def mat_ndim(self, value: int) -> None:
        """Set number of dimensions of the raw matrix"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def mat_shape(self) -> tuple:
        """Shape of the raw matrix."""
        raise NotImplementedError("property is abstract")
    @mat_shape.setter
    @abstractmethod
    def mat_shape(self, value: tuple) -> None:
        """Set the shape of the raw matrix"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    @abstractmethod
    def copy(
            self
        ) -> 'Matrix':
        """
        Make a shallow copy of the Matrix.

        Returns
        -------
        out : Matrix
            A shallow copy of the original Matrix.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def deepcopy(
            self, 
            memo: dict
        ) -> 'Matrix':
        """
        Make a deep copy of the Matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : Matrix
            A deep copy of the original Matrix.
        """
        raise NotImplementedError("method is abstract")

    ######### Matrix element copy-on-manipulation ##########
    @abstractmethod
    def adjoin(
            self, 
            values: Union['Matrix',numpy.ndarray], 
            axis: int, 
            **kwargs: dict
        ) -> 'Matrix':
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : Matrix or numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are adjoined.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    @abstractmethod
    def delete(
            self, 
            obj: Union[int,slice,Sequence], 
            axis: int, 
            **kwargs: dict
        ) -> 'Matrix':
        """
        Delete sub-arrays along an axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to delete the subarray defined by obj.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    @abstractmethod
    def insert(
            self, 
            obj: Union[int,slice,Sequence], 
            values: ArrayLike, 
            axis: int, 
            **kwargs: dict
        ) -> 'Matrix':
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : ArrayLike
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    @abstractmethod
    def select(
            self, 
            indices: ArrayLike, 
            axis: int, 
            **kwargs: dict
        ) -> 'Matrix':
        """
        Select certain values from the matrix.

        Parameters
        ----------
        indices : ArrayLike (Nj, ...)
            The indices of the values to select.
        axis : int
            The axis along which values are selected.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            The output matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    # TODO: there are discrepancies between this as a static method and other
    #       related methods such as concat_taxa being class methods. Consider
    #       converting this method to a class method as well.
    @classmethod
    @abstractmethod
    def concat(
            cls,
            mats: Sequence, 
            axis: int, 
            **kwargs: dict
        ) -> 'Matrix':
        """
        Concatenate matrices together along an axis.

        Parameters
        ----------
        mats : Sequence of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        axis : int
            The axis along which the arrays will be joined.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")



################################## Utilities ###################################
def check_is_Matrix(v: object, vname: str) -> None:
    """
    Check if object is of type Matrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, Matrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,Matrix.__name__,type(v).__name__))
