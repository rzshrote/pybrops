from . import Matrix

class MutableMatrix(Matrix):
    """docstring for MutableMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(MutableMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element in-place-manipulation #########
    def append(self, values, axis, **kwargs):
        """
        Append values to the Matrix.

        Parameters
        ----------
        values : numpy.ndarray
            Values are appended to append to the matrix.
        axis : int
            The axis along which values are appended.
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def remove(self, obj, axis, **kwargs):
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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_MutableMatrix(v):
    return isinstance(v, MutableMatrix)

def check_is_MutableMatrix(v, varname):
    if not isinstance(v, MutableMatrix):
        raise TypeError("'%s' must be a MutableMatrix." % varname)

def cond_check_is_MutableMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_MutableMatrix(v, varname)
