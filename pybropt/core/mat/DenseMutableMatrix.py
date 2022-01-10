from pybropt.core.mat.DenseMatrix import DenseMatrix
from pybropt.core.mat.DenseMatrix import is_DenseMatrix
from pybropt.core.mat.MutableMatrix import MutableMatrix

class DenseMutableMatrix(DenseMatrix,MutableMatrix):
    """docstring for DenseMutableMatrix."""

    def __init__(self, mat, **kwargs):
        """
        Constructor for DenseMutableMatrix

        Parameters
        ----------
        mat : numpy.ndarray
            Dense matrix to use to create object.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseMutableMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element in-place-manipulation #########
    def append(self, values, axis = -1, **kwargs):
        """
        Append values to the matrix.

        Parameters
        ----------
        values : DenseMatrix, numpy.ndarray
            Values are appended to append to the matrix.
        axis : int
            The axis along which values are appended.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a DenseMatrix extract Matrix.mat values
        if is_DenseMatrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseMatrix or numpy.ndarray")

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

    def incorp(self, obj, values, axis = -1, **kwargs):
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
        if is_DenseMatrix(values):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseMatrix or numpy.ndarray")

        # incorporate values
        self._mat = numpy.insert(self._mat, obj, values, axis)



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseMutableMatrix(v):
    """Test whether an object is a DenseMutableMatrix"""
    return isinstance(v, DenseMutableMatrix)

def check_is_DenseMutableMatrix(v, varname):
    """Raise error if object is not a DenseMutableMatrix"""
    if not isinstance(v, DenseMutableMatrix):
        raise TypeError("'%s' must be a DenseMutableMatrix." % varname)

def cond_check_is_DenseMutableMatrix(v, varname, cond=(lambda s: s is not None)):
    """Conditionally raise error if object is not a DenseMutableMatrix"""
    if cond(v):
        check_is_DenseMutableMatrix(v, varname)
