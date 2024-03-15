"""
Module implementing a dense matrix and associated error checking routines.
"""

__all__ = [
    "DenseMatrix",
    "check_is_DenseMatrix",
]

import copy
from pathlib import Path
import numpy
import h5py
from typing import Iterator, Optional, Sequence, Union
from numpy.typing import ArrayLike

from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group, check_h5py_File_is_readable, check_h5py_File_is_writable
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.util import get_axis
from pybrops.core.util.h5py import h5py_File_read_ndarray, h5py_File_write_dict

class DenseMatrix(Matrix):
    """
    A concrete class for dense matrices.
    Dense matrices utilize numpy.ndarray's for data storage.

    The purpose of this concrete class is to implement base functionality for:
        1) Dense matrix mathematical operators
        2) Dense matrix logical & bitwise operators
        3) Dense matrix container operators
        4) Dense matrix copy operators
        5) Dense matrix read-only matrix shape changing routines.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Construct a dense matrix object.

        Parameters
        ----------
        mat : numpy.ndarray
            Matrix to store.
        kwargs : dict
            Additional keyword arguments.
        """
        self.mat = mat

    ############## Forward numeric operators ###############
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
        return self._mat + value

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
        return self._mat - value

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
        return self._mat * value

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
        return self._mat @ value

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
        return self._mat / value

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
        return self._mat // value

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
        return self._mat % value

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
        return divmod(self._mat, value)

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
        return self._mat ** value

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
        return self._mat << value

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
        return self._mat >> value

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
        return self._mat & value

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
        return self._mat ^ value

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
        return self._mat | value

    ############# Backwards numeric operators ##############
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
        return value + self._mat

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
        return value - self._mat

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
        return value * self._mat

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
        return value @ self._mat

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
        return value / self._mat

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
        return value // self._mat

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
        return value % self._mat

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
        return divmod(value, self._mat)

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
        return value << self._mat

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
        return value >> self._mat

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
        return value & self._mat

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
        return value ^ self._mat

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
        return value | self._mat

    ############# Augmented numeric operators ##############
    def __iadd__(self, value: object) -> None:
        """
        Elementwise add assign matrices
        
        Parameters
        ----------
        value : object
            Object which to add.
        """
        self._mat += value

    def __isub__(self, value: object) -> None:
        """
        Elementwise subtract assign matrices
        
        Parameters
        ----------
        value : object
            Object which to subtract.
        """
        self._mat -= value

    def __imul__(self, value: object) -> None:
        """
        Elementwise multiply assign matrices
        
        Parameters
        ----------
        value : object
            Object which to multiply.
        """
        self._mat *= value

    def __imatmul__(self, value: object) -> None:
        """
        Multiply assign matrices
        
        Parameters
        ----------
        value : object
            Object which to matrix multiply.
        """
        self._mat @= value

    def __itruediv__(self, value: object) -> None:
        """
        Elementwise true divide assign matrices
        
        Parameters
        ----------
        value : object
            Object which to divide.
        """
        self._mat /= value

    def __ifloordiv__(self, value: object) -> None:
        """
        Elementwise floor divide assign matrices
        
        Parameters
        ----------
        value : object
            Object which to floor divide.
        """
        self._mat //= value

    def __imod__(self, value: object) -> None:
        """
        Elementwise modulus assign matrices
        
        Parameters
        ----------
        value : object
            Object which to modulo.
        """
        self._mat %= value

    def __ipow__(self, value: object) -> None:
        """
        Elementwise exponent assign matrices
        
        Parameters
        ----------
        value : object
            Object which to exponentiate.
        """
        self._mat **= value

    def __ilshift__(self, value: object) -> None:
        """
        Elementwise left shift assign matrices
        
        Parameters
        ----------
        value : object
            Object which to left shift.
        """
        self._mat <<= value

    def __irshift__(self, value: object) -> None:
        """
        Elementwise right shift assign matrices
        
        Parameters
        ----------
        value : object
            Object which to right shift.
        """
        self._mat >>= value

    def __iand__(self, value: object) -> None:
        """
        Elementwise bitwise AND assign matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise AND.
        """
        self._mat &= value

    def __ixor__(self, value: object) -> None:
        """
        Elementwise bitwise XOR assign matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise XOR.
        """
        self._mat ^= value

    def __ior__(self, value: object) -> None:
        """
        Elementwise bitwise OR assign matrices
        
        Parameters
        ----------
        value : object
            Object which to bitwise OR.
        """
        self._mat |= value

    ################## Logical operators ###################
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
        return self._mat < value

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
        return self._mat <= value

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
        return self._mat == value

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
        return self._mat != value

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
        return self._mat > value

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
        return self._mat >= value

    ################# Container operators ##################
    def __len__(self) -> int:
        """
        Get the length of the Matrix.

        Returns
        -------
        out : int
            The length of the Matrix.
        """
        return len(self._mat)

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
        return self._mat[key]

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
        self._mat[key] = value

    def __delitem__(self, key: object) -> None:
        """
        Delete an item in the Matrix.

        Parameters
        ----------
        key : object
            Key of the item which to delete.
        """
        del self._mat[key]

    def __iter__(self) -> Iterator:
        """
        Get an iterator for the Matrix.

        Returns
        -------
        out : Iterator
            An iterator for the Matrix.
        """
        return iter(self._mat)

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseMatrix
            A copy of the DenseMatrix.
        """
        return self.__class__(
            mat = copy.copy(self.mat)
        )

    def __deepcopy__(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseMatrix
            A deep copy of the DenseMatrix.
        """
        return self.__class__(
            mat = copy.deepcopy(self.mat, memo)
        )

    ########### Miscellaneous special functions ############
    def __repr__(
            self
        ) -> str:
        """
        Return repr(self).
        
        Returns
        -------
        out : str
            A representation of the object.
        """
        return "<{0} of shape {1} at {2}>".format(
            type(self).__name__,
            str(self.mat_shape),
            hex(id(self))
        )

    ############################ Object Properties #############################

    ##################### Matrix Data ######################
    @property
    def mat(self) -> numpy.ndarray:
        """Pointer to raw numpy.ndarray object."""
        return self._mat
    @mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set pointer to raw numpy.ndarray object."""
        # The only assumption is that mat is a numpy.ndarray matrix.
        # Let the user decide whether to overwrite error checks.
        check_is_ndarray(value, "mat")
        self._mat = value
    
    @property
    def mat_ndim(self) -> int:
        """Number of dimensions of the raw numpy.ndarray."""
        return self._mat.ndim
    @mat_ndim.setter
    def mat_ndim(self, value: int) -> None:
        """Set number of dimensions of the raw numpy.ndarray"""
        error_readonly("mat_ndim")
    
    @property
    def mat_shape(self) -> tuple:
        """Shape of the raw numpy.ndarray."""
        return self._mat.shape
    @mat_shape.setter
    def mat_shape(self, value: tuple) -> None:
        """Set the shape of the raw numpy.ndarray"""
        error_readonly("mat_shape")
    
    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseMatrix':
        """
        Make a shallow copy of the Matrix.

        Returns
        -------
        out : DenseMatrix
            A shallow copy of the original DenseMatrix.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseMatrix':
        """
        Make a deep copy of the Matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseMatrix
            A deep copy of the original DenseMatrix.
        """
        return copy.deepcopy(self, memo)

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(
            self, 
            values: Union['DenseMatrix',numpy.ndarray], 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseMatrix':
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : DenseMatrix, numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are appended.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseMatrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a Matrix extract Matrix.mat values
        if isinstance(values, DenseMatrix):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseMatrix or numpy.ndarray")

        # append values
        mat = numpy.append(self._mat, values, axis)

        # create new output
        out = self.__class__(mat = mat, **kwargs)

        return out

    def delete(
            self, 
            obj: Union[int,slice,Sequence], 
            axis: int, 
            **kwargs: dict
        ) -> 'DenseMatrix':
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
        out : DenseMatrix
            A ``DenseMatrix`` with deleted elements. Note that concat does not occur
            in-place: a new ``DenseMatrix`` is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # delete values
        mat = numpy.delete(self._mat, obj, axis)

        # create new output
        out = self.__class__(mat = mat, **kwargs)

        return out

    def insert(
            self, 
            obj: Union[int,slice,Sequence], 
            values: ArrayLike, 
            axis: int, 
            **kwargs: dict
        ) -> 'DenseMatrix':
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
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
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a Matrix extract Matrix.mat values
        if isinstance(values, DenseMatrix):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseMatrix or numpy.ndarray")

        # append values
        mat = numpy.insert(self._mat, obj, values, axis)

        # create new output
        out = self.__class__(mat = mat, **kwargs)

        return out

    def select(
            self, 
            indices: ArrayLike, 
            axis: int, 
            **kwargs: dict
        ) -> 'DenseMatrix':
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
        out : DenseMatrix
            A ``DenseMatrix`` with values selected. Note that select does not
            occur in-place: a new ``DenseMatrix`` is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # select values
        mat = numpy.take(self._mat, indices, axis)

        # create new output
        out = self.__class__(mat = mat, **kwargs)

        return out

    @staticmethod
    def concat(
        mats: ArrayLike, 
        axis: int, 
        **kwargs: dict
    ) -> 'DenseMatrix':
        """
        Concatenate matrices together along an axis.

        Parameters
        ----------
        mats : ArrayLike of DenseMatrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        axis : int
            The axis along which the arrays will be joined.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            A concatenated ``DenseMatrix``. Note that concat does not occur in-place:
            a new ``DenseMatrix`` is allocated and filled.
        """
        # gather raw matrices
        mat_tp = tuple(m.mat for m in mats)

        # concatenate matrices along axis
        mat = numpy.concatenate(mat_tp, axis)

        # create new output
        out = mats[0].__class__(mat = mat, **kwargs)

        return out

    ################### Matrix File I/O ####################
    def to_hdf5(
            self, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            overwrite: bool = True,
        ) -> None:
        """
        Write GenotypeMatrix to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name to which to write. File is closed after writing.
            If ``h5py.File``, an opened HDF5 file to which to write. File is not closed after writing.

        groupname : str, None
            If ``str``, an HDF5 group name under which the ``DenseMatrix`` data is stored.
            If ``None``, the ``DenseMatrix`` is written to the base HDF5 group.

        overwrite : bool
            Whether to overwrite values in an HDF5 file if a field already exists.
        """
        ########################################################
        ############ process ``filename`` argument #############

        # HDF5 file object
        h5file = None

        # if we have a string or Path, open HDF5 file in append (``r+``) mode
        if isinstance(filename, (str,Path)):
            h5file = h5py.File(filename, "a")

        # elif we have an h5py.File, make sure mode is writable, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_writable(filename)
            h5file = filename
        
        # else raise TypeError
        else:
            raise TypeError(
                "``filename`` must be of type ``str``, ``Path``, or ``h5py.File`` but received type ``{0}``".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        #### write data to HDF5 file and (optionally) close ####

        # data dictionary
        data = {
            "mat" : self.mat,
        }
        
        # save data to file
        h5py_File_write_dict(h5file, groupname, data, overwrite)
        
        # close the file, only if the provided filename was a string and not a h5py.File.
        if isinstance(filename, str):
            h5file.close()

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None
        ) -> 'DenseMatrix':
        """
        Read DenseMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which DenseMatrix data is stored.
            If None, DenseMatrix is read from base HDF5 group.

        Returns
        -------
        out : DenseMatrix
            A dense matrix read from file.
        """
        ########################################################
        ############ process ``filename`` argument #############

        # HDF5 file object
        h5file = None

        # if we have a string or Path, open HDF5 file in read (``r``) mode
        if isinstance(filename, (str,Path)):
            check_file_exists(filename)
            h5file = h5py.File(filename, "r")

        # elif we have an ``h5py.File``, make sure mode is in at least ``r`` mode, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_readable(filename)
            h5file = filename
        
        # else raise TypeError
        else:
            raise TypeError(
                "``filename`` must be of type ``str``, ``Path``, or ``h5py.File`` but received type ``{0}``".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # FIXME: errors if groupname == "" or "/"
            # if the group does not exist in the file, close and raise error
            check_h5py_File_has_group(h5file, groupname)

            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        ######## check that we have all required fields ########

        # all required arguments
        required_fields = ["mat"]

        # for each required field, check if the field exists in the HDF5 file.
        for field in required_fields:
            check_h5py_File_has_group(h5file, groupname + field)

        ########################################################
        ### read data from HDF5 file and (optionally) close ####
        
        # output dictionary
        data = {
            "mat": None
        }
        
        ##################################
        ### read mandatory data fields ###

        # read mat array (ndarray dtype = int8)
        data["mat"] = h5py_File_read_ndarray(h5file, groupname + "mat")
        
        ######################
        ### close the file ###

        # close the file, only if the provided fieldname was a string an not an h5py.File.
        if isinstance(filename, str):
            h5file.close()
        
        ########################################################
        ################### Object creation ####################
        
        # create object from read data
        out = cls(
            mat = data["mat"],
        )

        return out



################################## Utilities ###################################
def check_is_DenseMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,DenseMatrix.__name__,type(v).__name__))
