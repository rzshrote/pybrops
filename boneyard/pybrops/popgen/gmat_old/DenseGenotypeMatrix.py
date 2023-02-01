import copy
import numpy

from . import GenotypeMatrix
from pybrops.core.error import check_is_ndarray
from pybrops.core.mat import DenseMutableMatrix
from pybrops.core.mat import get_axis
from pybrops.core.mat import is_Matrix

class DenseGenotypeMatrix(DenseMutableMatrix,GenotypeMatrix):
    """docstring for DenseGenotypeMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, **kwargs: dict):
        super(DenseGenotypeMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

    #################### Matrix copying ####################
    def __copy__(self):
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : Matrix
        """
        return self.__class__(
            mat = copy.copy(self.mat)
        )

    def __deepcopy__(self, memo):
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict

        Returns
        -------
        out : Matrix
        """
        return self.__class__(
            mat = copy.deepcopy(self.mat)
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    #################### Genotypic Data ####################
    def mat():
        doc = "The mat property."
        def fget(self):
            return self._mat
        def fset(self, value):
            # The only assumption is that mat is a numpy.ndarray matrix.
            # Let the user decide whether to overwrite error checks.
            check_is_ndarray(value, "mat")
            self._mat = value
        def fdel(self):
            del self._mat
        return locals()
    mat = property(**mat())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    #################### Matrix copying ####################
    def copy(self):
        """
        Make a shallow copy of the Matrix.

        Returns
        -------
        out : Matrix
            A shallow copy of the original Matrix.
        """
        return copy.copy(self)

    def deepcopy(self, memo = None):
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
        return copy.deepcopy(self, memo)

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(self, values, axis = -1, **kwargs: dict):
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : Matrix or numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are appended.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a Matrix extract Matrix.mat values
        if is_Matrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type Matrix or numpy.ndarray")

        # append values
        mat = numpy.append(self._mat, values, axis)

        # create new output
        out = DenseGenotypeMatrix(mat = mat)

        return out

    def delete(self, obj, axis = -1, **kwargs: dict):
        """
        Delete sub-arrays along an axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to delete the subarray defined by obj.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # delete values
        mat = numpy.delete(self._mat, obj, axis)

        # create new output
        out = DenseGenotypeMatrix(mat = mat)

        return out

    def insert(self, obj, values, axis = -1, **kwargs: dict):
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        **kwargs
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
        if is_Matrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type Matrix or numpy.ndarray")

        # append values
        mat = numpy.insert(self._mat, obj, values, axis)

        # create new output
        out = DenseGenotypeMatrix(mat = mat)

        return out

    def select(self, indices, axis = -1, **kwargs: dict):
        """
        Select certain values from the matrix.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        axis : int
            The axis along which values are selected.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            The output matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # select values
        mat = numpy.take(self._mat, indices, axis)

        # create new output
        out = DenseGenotypeMatrix(mat = mat)

        return out

    @staticmethod
    def concat(mats, axis = -1, **kwargs: dict):
        """
        Concatenate matrices together along an axis.

        Parameters
        ----------
        mats : array_like of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        axis : int
            The axis along which the arrays will be joined.
        **kwargs
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        # gather raw matrices
        rawmat = tuple(m.mat for m in mats)

        # concatenate matrices along axis
        mat = numpy.concatenate(rawmat, axis)

        # create new output
        out = DenseGenotypeMatrix(mat = mat)

        return out

    ######### Matrix element in-place-manipulation #########
    def append(self, values, axis = -1, **kwargs: dict):
        """
        Append values to the matrix.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        axis : int
            The axis along which values are appended.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a Matrix extract Matrix.mat values
        if is_Matrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type Matrix or numpy.ndarray")

        # append values
        self._mat = numpy.append(self._mat, values, axis)

    def remove(self, obj, axis = -1, **kwargs: dict):
        """
        Remove sub-arrays along an axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to remove the subarray defined by obj.
        **kwargs
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # delete values
        self._mat = numpy.delete(self._mat, obj, axis)

    def incorp(self, obj, values, axis = -1, **kwargs: dict):
        """
        Incorporate values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : array_like
            Values to incorporate into the matrix.
        axis : int
            The axis along which values are incorporated.
        **kwargs
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a Matrix extract Matrix.mat values
        if is_Matrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseMatrix or numpy.ndarray")

        # incorporate values
        self._mat = numpy.insert(self._mat, obj, values, axis)



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseGenotypeMatrix(v):
    return isinstance(v, DenseGenotypeMatrix)

def check_is_DenseGenotypeMatrix(v, varname):
    if not isinstance(v, DenseGenotypeMatrix):
        raise TypeError("'%s' must be a DenseGenotypeMatrix." % varname)

def cond_check_is_DenseGenotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DenseGenotypeMatrix(v, varname)
