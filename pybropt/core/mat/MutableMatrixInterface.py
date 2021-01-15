class MutableMatrixInterface:
    """docstring for MutableMatrixInterface."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(MutableMatrixInterface, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############# Matrix element manipulation ##############
    def append(self, values, axis, **kwargs):
        """
        Append values to the matrix.

        Parameters
        ----------
        values : numpy.ndarray
            Values are appended to append to the matrix.
        axis : int
            The axis along which values are appended.
        """
        raise NotImplementedError("method is abstract")

    def delete(self, obj, axis, **kwargs):
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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")


################################################################################
################################## Utilities ###################################
################################################################################
def is_MutableMatrixInterface(v):
    return isinstance(v, MutableMatrixInterface)

def check_is_MutableMatrixInterface(v, varname):
    if not isinstance(v, MutableMatrixInterface):
        raise TypeError("'%s' must be a MutableMatrixInterface." % varname)

def cond_check_is_MutableMatrixInterface(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_MutableMatrixInterface(v, varname)
