class Matrix:
    """An abstract class for matrix wrapper objects."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Matrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(Matrix, self).__init__()

    ############## Forward numeric operators ###############
    def __add__(self, value):
        raise NotImplementedError("method is abstract")

    def __sub__(self, value):
        raise NotImplementedError("method is abstract")

    def __mul__(self, value):
        raise NotImplementedError("method is abstract")

    def __matmul__(self, value):
        raise NotImplementedError("method is abstract")

    def __truediv__(self, value):
        raise NotImplementedError("method is abstract")

    def __floordiv__(self, value):
        raise NotImplementedError("method is abstract")

    def __mod__(self, value):
        raise NotImplementedError("method is abstract")

    def __divmod__(self, value):
        raise NotImplementedError("method is abstract")

    def __pow__(self, value):
        raise NotImplementedError("method is abstract")

    def __lshift__(self, value):
        raise NotImplementedError("method is abstract")

    def __rshift__(self, value):
        raise NotImplementedError("method is abstract")

    def __and__(self, value):
        raise NotImplementedError("method is abstract")

    def __xor__(self, value):
        raise NotImplementedError("method is abstract")

    def __or__(self, value):
        raise NotImplementedError("method is abstract")

    ############# Backwards numeric operators ##############
    def __radd__(self, value):
        raise NotImplementedError("method is abstract")

    def __rsub__(self, value):
        raise NotImplementedError("method is abstract")

    def __rmul__(self, value):
        raise NotImplementedError("method is abstract")

    def __rmatmul__(self, value):
        raise NotImplementedError("method is abstract")

    def __rtruediv__(self, value):
        raise NotImplementedError("method is abstract")

    def __rfloordiv__(self, value):
        raise NotImplementedError("method is abstract")

    def __rmod__(self, value):
        raise NotImplementedError("method is abstract")

    def __rdivmod__(self, value):
        raise NotImplementedError("method is abstract")

    def __rlshift__(self, value):
        raise NotImplementedError("method is abstract")

    def __rrshift__(self, value):
        raise NotImplementedError("method is abstract")

    def __rand__(self, value):
        raise NotImplementedError("method is abstract")

    def __rxor__(self, value):
        raise NotImplementedError("method is abstract")

    def __ror__(self, value):
        raise NotImplementedError("method is abstract")

    ############# Augmented numeric operators ##############
    def __iadd__(self, value):
        raise NotImplementedError("method is abstract")

    def __isub__(self, value):
        raise NotImplementedError("method is abstract")

    def __imul__(self, value):
        raise NotImplementedError("method is abstract")

    def __imatmul__(self, value):
        raise NotImplementedError("method is abstract")

    def __itruediv__(self, value):
        raise NotImplementedError("method is abstract")

    def __ifloordiv__(self, value):
        raise NotImplementedError("method is abstract")

    def __imod__(self, value):
        raise NotImplementedError("method is abstract")

    def __ipow__(self, value):
        raise NotImplementedError("method is abstract")

    def __ilshift__(self, value):
        raise NotImplementedError("method is abstract")

    def __irshift__(self, value):
        raise NotImplementedError("method is abstract")

    def __iand__(self, value):
        raise NotImplementedError("method is abstract")

    def __ixor__(self, value):
        raise NotImplementedError("method is abstract")

    def __ior__(self, value):
        raise NotImplementedError("method is abstract")

    ################## Logical operators ###################
    def __lt__(self, value):
        raise NotImplementedError("method is abstract")

    def __le__(self, value):
        raise NotImplementedError("method is abstract")

    def __eq__(self, value):
        raise NotImplementedError("method is abstract")

    def __ne__(self, value):
        raise NotImplementedError("method is abstract")

    def __gt__(self, value):
        raise NotImplementedError("method is abstract")

    def __ge__(self, value):
        raise NotImplementedError("method is abstract")

    ################# Container operators ##################
    def __len__(self):
        raise NotImplementedError("method is abstract")

    def __getitem__(self, key):
        """
        Get a specific key within the raw underlying matrix.
        """
        raise NotImplementedError("method is abstract")

    def __setitem__(self, key, value):
        """
        Set a specific key within the raw underlying matrix.
        """
        raise NotImplementedError("method is abstract")

    def __delitem__(self, key):
        """
        Delete a specific key within the raw underlying matrix.
        """
        raise NotImplementedError("method is abstract")

    def __iter__(self):
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def mat():
        doc = "Pointer to raw matrix object."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    mat = property(**mat())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############# Matrix element manipulation ##############
    def append(self, values, axis, **kwargs):
        """
        Append values to the matrix.

        Parameters
        ----------
        values : array_like
            Values are appended to append to the matrix.
        axis : int
            The axis along which values are appended.
        **kwargs
            Additional keyword arguments.
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



################################################################################
################################## Utilities ###################################
################################################################################
def is_Matrix(v):
    return isinstance(v, Matrix)

def check_is_Matrix(v, varname):
    if not isinstance(v, Matrix):
        raise TypeError("'%s' must be a Matrix." % varname)

def cond_check_is_Matrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_Matrix(v, varname)
