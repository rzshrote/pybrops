import copy
import numpy

from . import GenotypeMatrix
from pybropt.core.error import check_is_ndarray
from pybropt.core.mat import get_axis

class DenseGenotypeMatrix(GenotypeMatrix):
    """docstring for DenseGenotypeMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, **kwargs):
        super(DenseGenotypeMatrix, self).__init__(
            **kwargs
        )
        self.mat = mat

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
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : Matrix
        """
        return self.__copy__()

    ############# Matrix element manipulation ##############
    def append(self, values, axis = -1, **kwargs):
        """
        Append values to the matrix.

        Parameters
        ----------
        values : numpy.ndarray
            Values are appended to append to the matrix.
        axis : int
            The axis along which values are appended.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)
        # append values
        self._mat = numpy.append(self._mat, values, axis)

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
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)
        # delete values
        self._mat = numpy.delete(self._mat, obj, axis)

    def insert(self, obj, values, axis, **kwargs):
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : array_like
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        **kwargs
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)
        # insert values
        self._mat = numpy.insert(self._mat, obj, values, axis)

    def select(self, obj, axis, **kwargs):
        """
        Select certain values from the GenotypeMatrix.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices where values are selected.
        axis : int
            The axis along which values are selected.
        **kwargs
            Additional keyword arguments.
        """
        ndim = self._mat.ndim           # get number of dimensions
        axis = get_axis(axis, ndim)     # get axis

        # construct selection tuple
        sel = tuple(slice(None) if e != axis else obj for e in range(ndim))

        # select values
        smat = self._mat[sel]

        # create selection matrix
        sdgmat = DenseGenotypeMatrix(
            mat = mat,
        )

        return sdgmat



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
