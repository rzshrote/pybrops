import copy
import numpy
from pybrops.core.mat.Matrix import Matrix

from pybrops.core.error import error_readonly
from pybrops.core.error import check_is_ndarray

class DenseMatrix(Matrix):
    """
    DenseMatrix class.

    The purpose of this class is to provide dense matrix math operators, and
    copy on manipulation subroutines.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, **kwargs):
        """
        Parameters
        ----------
        mat : numpy.ndarray
            Matrix to store
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseMatrix, self).__init__(**kwargs)
        self.mat = mat

    ############## Forward numeric operators ###############
    def __add__(self, value):
        """Elementwise add matrices"""
        return self._mat + value

    def __sub__(self, value):
        """Elementwise subtract matrices"""
        return self._mat - value

    def __mul__(self, value):
        """Elementwise multiply matrices"""
        return self._mat * value

    def __matmul__(self, value):
        """Multiply matrices"""
        return self._mat @ value

    def __truediv__(self, value):
        """Elementwise divide matrices"""
        return self._mat / value

    def __floordiv__(self, value):
        """Elementwise floor divide matrices"""
        return self._mat // value

    def __mod__(self, value):
        """Elementwise modulus matrices"""
        return self._mat % value

    def __divmod__(self, value):
        """Elementwise divmod matrices"""
        return divmod(self._mat, value)

    def __pow__(self, value):
        """Elementwise exponent matrices"""
        return self._mat ** value

    def __lshift__(self, value):
        """Elementwise bitwise left shift matrices"""
        return self._mat << value

    def __rshift__(self, value):
        """Elementwise bitwise right shift matrices"""
        return self._mat >> value

    def __and__(self, value):
        """Elementwise bitwise and matrices"""
        return self._mat & value

    def __xor__(self, value):
        """Elementwise bitwise xor matrices"""
        return self._mat ^ value

    def __or__(self, value):
        """Elementwise bitwise or matrices"""
        return self._mat | value

    ############# Backwards numeric operators ##############
    def __radd__(self, value):
        """Reverse elementwise add matrices"""
        return value + self._mat

    def __rsub__(self, value):
        """Reverse elementwise subtract matrices"""
        return value - self._mat

    def __rmul__(self, value):
        """Reverse elementwise multiply matrices"""
        return value * self._mat

    def __rmatmul__(self, value):
        """Reverse multiply matrices"""
        return value @ self._mat

    def __rtruediv__(self, value):
        """Reverse elementwise divide matrices"""
        return value / self._mat

    def __rfloordiv__(self, value):
        """Reverse elementwise floor divide matrices"""
        return value // self._mat

    def __rmod__(self, value):
        """Reverse elementwise modulus matrices"""
        return value % self._mat

    def __rdivmod__(self, value):
        """Reverse elementwise divmod matrices"""
        return divmod(value, self._mat)

    def __rlshift__(self, value):
        """Reverse elementwise bitwise left shift matrices"""
        return value << self._mat

    def __rrshift__(self, value):
        """Reverse elementwise bitwise right shift matrices"""
        return value >> self._mat

    def __rand__(self, value):
        """Reverse elementwise bitwise and matrices"""
        return value & self._mat

    def __rxor__(self, value):
        """Reverse elementwise bitwise xor matrices"""
        return value ^ self._mat

    def __ror__(self, value):
        """Reverse elementwise bitwise or matrices"""
        return value | self._mat

    ############# Augmented numeric operators ##############
    def __iadd__(self, value):
        """Elementwise add assign matrices"""
        self._mat += value

    def __isub__(self, value):
        """Elementwise subtract assign matrices"""
        self._mat -= value

    def __imul__(self, value):
        """Elementwise multiply assign matrices"""
        self._mat *= value

    def __imatmul__(self, value):
        """Multiply assign matrices"""
        self._mat @= value

    def __itruediv__(self, value):
        """Elementwise true divide assign matrices"""
        self._mat /= value

    def __ifloordiv__(self, value):
        """Elementwise floor divide assign matrices"""
        self._mat //= value

    def __imod__(self, value):
        """Elementwise modulus assign matrices"""
        self._mat %= value

    def __ipow__(self, value):
        """Elementwise exponent assign matrices"""
        self._mat **= value

    def __ilshift__(self, value):
        """Elementwise left shift assign matrices"""
        self._mat <<= value

    def __irshift__(self, value):
        """Elementwise right shift assign matrices"""
        self._mat >>= value

    def __iand__(self, value):
        """Elementwise bitwise and assign matrices"""
        self._mat &= value

    def __ixor__(self, value):
        """Elementwise bitwise xor assign matrices"""
        self._mat ^= value

    def __ior__(self, value):
        """Elementwise bitwise or assign matrices"""
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
            mat = copy.deepcopy(self.mat, memo)
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ##################### Matrix Data ######################
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

    def mat_ndim():
        doc = "Number of dimensions of the raw matrix property."
        def fget(self):
            """Get number of dimensions of the raw matrix"""
            return self._mat.ndim
        def fset(self, value):
            """Set number of dimensions of the raw matrix"""
            error_readonly("mat_ndim")
        def fdel(self):
            """Delete number of dimensions of the raw matrix"""
            error_readonly("mat_ndim")
        return locals()
    mat_ndim = property(**mat_ndim())

    def mat_shape():
        doc = "Shape of the raw matrix property."
        def fget(self):
            """Get the shape of the raw matrix"""
            return self._mat.shape
        def fset(self, value):
            """Set the shape of the raw matrix"""
            error_readonly("mat_shape")
        def fdel(self):
            """Delete the shape of the raw matrix"""
            error_readonly("mat_shape")
        return locals()
    mat_shape = property(**mat_shape())

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
    def adjoin(self, values, axis = -1, **kwargs):
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : DenseMatrix or numpy.ndarray
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
        if is_DenseMatrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseMatrix or numpy.ndarray")

        # append values
        mat = numpy.append(self._mat, values, axis)

        # create new output
        out = self.__class__(mat = mat, **kwargs)

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
        kwargs : dict
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
        out = self.__class__(mat = mat, **kwargs)

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
        if is_DenseMatrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseMatrix or numpy.ndarray")

        # append values
        mat = numpy.insert(self._mat, obj, values, axis)

        # create new output
        out = self.__class__(mat = mat, **kwargs)

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
        kwargs : dict
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
        out = self.__class__(mat = mat, **kwargs)

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
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        # gather raw matrices
        mat_tp = tuple(m.mat for m in mats)

        # concatenate matrices along axis
        mat = numpy.concatenate(mat_tp, axis)

        # create new output
        out = self.__class__(mat = mat, **kwargs)

        return out



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseMatrix(v):
    """Return whether an object is a DenseMatrix or not"""
    return isinstance(v, DenseMatrix)

def check_is_DenseMatrix(v, varname):
    """Raise TypeError if object is not a DenseMatrix"""
    if not isinstance(v, DenseMatrix):
        raise TypeError("'%s' must be a DenseMatrix." % varname)

def cond_check_is_DenseMatrix(v, varname, cond=(lambda s: s is not None)):
    """If object is not None, raise TypeError if object is not a DenseMatrix"""
    if cond(v):
        check_is_DenseMatrix(v, varname)