import copy
import numpy

from . import BreedingValueMatrix
from pybropt.core.mat import get_axis
from pybropt.core.mat import is_Matrix
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import cond_check_is_ndarray

class DenseBreedingValueMatrix(BreedingValueMatrix):
    """
    Partial implementation of the BreedingValueMatrix interface.
    Implements matrix numeric, logical, and container operators.
    Implements matrix checking.
    All else is abstract.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, **kwargs):
        """
        BreedingValueMatrix constructor

        Parameters
        ----------
        mat : numpy.ndarray
            A float64 matrix of breeding values of shape (n, t).
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.

        """
        super(DenseBreedingValueMatrix, self).__init__(**kwargs)
        self.mat = mat

    ############## Forward numeric operators ###############
    def __add__(self, value):
        return self._mat + value

    def __sub__(self, value):
        return self._mat - value

    def __mul__(self, value):
        return self._mat * value

    def __matmul__(self, value):
        return self._mat @ value

    def __truediv__(self, value):
        return self._mat / value

    def __floordiv__(self, value):
        return self._mat // value

    def __mod__(self, value):
        return self._mat % value

    def __divmod__(self, value):
        return divmod(self._mat, value)

    def __pow__(self, value):
        return self._mat ** value

    def __lshift__(self, value):
        return self._mat << value

    def __rshift__(self, value):
        return self._mat >> value

    def __and__(self, value):
        return self._mat & value

    def __xor__(self, value):
        return self._mat ^ value

    def __or__(self, value):
        return self._mat | value

    ############# Backwards numeric operators ##############
    def __radd__(self, value):
        return value + self._mat

    def __rsub__(self, value):
        return value - self._mat

    def __rmul__(self, value):
        return value * self._mat

    def __rmatmul__(self, value):
        return value @ self._mat

    def __rtruediv__(self, value):
        return value / self._mat

    def __rfloordiv__(self, value):
        return value // self._mat

    def __rmod__(self, value):
        return value % self._mat

    def __rdivmod__(self, value):
        return divmod(value, self._mat)

    def __rlshift__(self, value):
        return value << self._mat

    def __rrshift__(self, value):
        return value >> self._mat

    def __rand__(self, value):
        return value & self._mat

    def __rxor__(self, value):
        return value ^ self._mat

    def __ror__(self, value):
        return value | self._mat

    ############# Augmented numeric operators ##############
    def __iadd__(self, value):
        self._mat += value

    def __isub__(self, value):
        self._mat -= value

    def __imul__(self, value):
        self._mat *= value

    def __imatmul__(self, value):
        self._mat @= value

    def __itruediv__(self, value):
        self._mat /= value

    def __ifloordiv__(self, value):
        self._mat //= value

    def __imod__(self, value):
        self._mat %= value

    def __ipow__(self, value):
        self._mat **= value

    def __ilshift__(self, value):
        self._mat <<= value

    def __irshift__(self, value):
        self._mat >>= value

    def __iand__(self, value):
        self._mat &= value

    def __ixor__(self, value):
        self._mat ^= value

    def __ior__(self, value):
        self._mat |= value

    ################## Logical operators ###################
    def __lt__(self, value):
        return self._mat < value

    def __le__(self, value):
        return self._mat <= value

    def __eq__(self, value):
        return self._mat == value

    def __ne__(self, value):
        return self._mat != value

    def __gt__(self, value):
        return self._mat > value

    def __ge__(self, value):
        return self._mat >= value

    ################# Container operators ##################
    def __len__(self):
        return len(self._mat)

    def __getitem__(self, key):
        return self._mat[key]

    def __setitem__(self, key, value):
        self._mat[key] = value

    def __delitem__(self, key):
        del self._mat[key]

    def __iter__(self):
        return iter(self._mat)

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

    ################# Breeding Value Data ##################
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
        return self.__copy__()

    def deepcopy(self, memo):
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
        return self.__deepcopy__(memo)

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(self, values, axis = -1, **kwargs):
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
        out = self.__class__(mat = mat)

        return out

    def delete(self, obj, axis = -1, **kwargs):
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
        out = self.__class__(mat = mat)

        return out

    def insert(self, obj, values, axis = -1, **kwargs):
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
        out = self.__class__(mat = mat)

        return out

    def select(self, indices, axis = -1, **kwargs):
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
        out = self.__class__(mat = mat)

        return out

    @staticmethod
    def concat(mats, axis = -1, **kwargs):
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
    def append(self, values, axis = -1, **kwargs):
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

    def remove(self, obj, axis = -1, **kwargs):
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

    def incorp(self, obj, values, axis, **kwargs):
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

        # incorporate values
        self._mat = numpy.insert(self._mat, obj, values, axis)

    # TODO: append, delete, insert, select


################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseBreedingValueMatrix(v):
    return isinstance(v, DenseBreedingValueMatrix)

def check_is_DenseBreedingValueMatrix(v, varname):
    if not isinstance(v, DenseBreedingValueMatrix):
        raise TypeError("'%s' must be a DenseBreedingValueMatrix." % varname)

def cond_check_is_DenseBreedingValueMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DenseBreedingValueMatrix(v, varname)
