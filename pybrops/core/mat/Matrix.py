"""
Module defining basal Matrix interfaces and associated error checking routines.
"""

class Matrix:
    """
    An abstract class for matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix mathematical operators
        2) Matrix logical & bitwise operators
        3) Matrix container operators
        4) Matrix copy operators
        5) Matrix read-only matrix shape changing routines.

    The shape of a Matrix should be immutable.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Matrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(Matrix, self).__init__()

    ############## Forward numeric operators ###############
    def __add__(self, value):
        """
        Return self + value.
        """
        raise NotImplementedError("method is abstract")

    def __sub__(self, value):
        """
        Return self - value.
        """
        raise NotImplementedError("method is abstract")

    def __mul__(self, value):
        """
        Return self * value.
        """
        raise NotImplementedError("method is abstract")

    def __matmul__(self, value):
        """
        Return self @ value.
        """
        raise NotImplementedError("method is abstract")

    def __truediv__(self, value):
        """
        Return self / value.
        """
        raise NotImplementedError("method is abstract")

    def __floordiv__(self, value):
        """
        Return self // value.
        """
        raise NotImplementedError("method is abstract")

    def __mod__(self, value):
        """
        Return self % value.
        """
        raise NotImplementedError("method is abstract")

    def __divmod__(self, value):
        """
        Return divmod(self, value).
        """
        raise NotImplementedError("method is abstract")

    def __pow__(self, value):
        """
        Return self ** value.
        """
        raise NotImplementedError("method is abstract")

    def __lshift__(self, value):
        """
        Return self << value.
        """
        raise NotImplementedError("method is abstract")

    def __rshift__(self, value):
        """
        Return self >> value.
        """
        raise NotImplementedError("method is abstract")

    def __and__(self, value):
        """
        Return self & value.
        """
        raise NotImplementedError("method is abstract")

    def __xor__(self, value):
        """
        Return self ^ value.
        """
        raise NotImplementedError("method is abstract")

    def __or__(self, value):
        """
        Return self | value.
        """
        raise NotImplementedError("method is abstract")

    ############# Backwards numeric operators ##############
    def __radd__(self, value):
        """
        Return value + self.
        """
        raise NotImplementedError("method is abstract")

    def __rsub__(self, value):
        """
        Return value - self.
        """
        raise NotImplementedError("method is abstract")

    def __rmul__(self, value):
        """
        Return value * self.
        """
        raise NotImplementedError("method is abstract")

    def __rmatmul__(self, value):
        """
        Return value @ self.
        """
        raise NotImplementedError("method is abstract")

    def __rtruediv__(self, value):
        """
        Return value / self.
        """
        raise NotImplementedError("method is abstract")

    def __rfloordiv__(self, value):
        """
        Return value // self.
        """
        raise NotImplementedError("method is abstract")

    def __rmod__(self, value):
        """
        Return value % self.
        """
        raise NotImplementedError("method is abstract")

    def __rdivmod__(self, value):
        """
        Return divmod(value, self).
        """
        raise NotImplementedError("method is abstract")

    def __rlshift__(self, value):
        """
        Return value << self.
        """
        raise NotImplementedError("method is abstract")

    def __rrshift__(self, value):
        """
        Return value >> self.
        """
        raise NotImplementedError("method is abstract")

    def __rand__(self, value):
        """
        Return value & self.
        """
        raise NotImplementedError("method is abstract")

    def __rxor__(self, value):
        """
        Return value ^ self.
        """
        raise NotImplementedError("method is abstract")

    def __ror__(self, value):
        """
        Return value | self.
        """
        raise NotImplementedError("method is abstract")

    ############# Augmented numeric operators ##############
    def __iadd__(self, value):
        """
        Return self += value.
        """
        raise NotImplementedError("method is abstract")

    def __isub__(self, value):
        """
        Return self -= value.
        """
        raise NotImplementedError("method is abstract")

    def __imul__(self, value):
        """
        Return self *= value.
        """
        raise NotImplementedError("method is abstract")

    def __imatmul__(self, value):
        """
        Return self @= value.
        """
        raise NotImplementedError("method is abstract")

    def __itruediv__(self, value):
        """
        Return self /= value.
        """
        raise NotImplementedError("method is abstract")

    def __ifloordiv__(self, value):
        """
        Return self //= value.
        """
        raise NotImplementedError("method is abstract")

    def __imod__(self, value):
        """
        Return self %= value.
        """
        raise NotImplementedError("method is abstract")

    def __ipow__(self, value):
        """
        Return self **= value.
        """
        raise NotImplementedError("method is abstract")

    def __ilshift__(self, value):
        """
        Return self <<= value.
        """
        raise NotImplementedError("method is abstract")

    def __irshift__(self, value):
        """
        Return self >>= value.
        """
        raise NotImplementedError("method is abstract")

    def __iand__(self, value):
        """
        Return self &= value.
        """
        raise NotImplementedError("method is abstract")

    def __ixor__(self, value):
        """
        Return self ^= value.
        """
        raise NotImplementedError("method is abstract")

    def __ior__(self, value):
        """
        Return self |= value.
        """
        raise NotImplementedError("method is abstract")

    ################## Logical operators ###################
    def __lt__(self, value):
        """
        Return self < value.
        """
        raise NotImplementedError("method is abstract")

    def __le__(self, value):
        """
        Return self <= value.
        """
        raise NotImplementedError("method is abstract")

    def __eq__(self, value):
        """
        Return self == value.
        """
        raise NotImplementedError("method is abstract")

    def __ne__(self, value):
        """
        Return self != value.
        """
        raise NotImplementedError("method is abstract")

    def __gt__(self, value):
        """
        Return self > value.
        """
        raise NotImplementedError("method is abstract")

    def __ge__(self, value):
        """
        Return self >= value.
        """
        raise NotImplementedError("method is abstract")

    ################# Container operators ##################
    def __len__(self):
        """
        Get the length of the raw underlying matrix.
        """
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
        """
        Create an iterator of the raw underlying matrix.
        """
        raise NotImplementedError("method is abstract")

    #################### Matrix copying ####################
    def __copy__(self):
        """
        Make a shallow copy of the Matrix.

        Returns
        -------
        out : Matrix
            A shallow copy of the original Matrix.
        """
        raise NotImplementedError("method is abstract")

    def __deepcopy__(self, memo):
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

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def mat():
        doc = "Pointer to raw matrix object."
        def fget(self):
            """Get pointer to raw matrix object"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set pointer to raw matrix object"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete raw matrix object"""
            raise NotImplementedError("method is abstract")
        return locals()
    mat = property(**mat())

    def mat_ndim():
        doc = "Number of dimensions of the raw matrix property."
        def fget(self):
            """Get number of dimensions of the raw matrix"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set number of dimensions of the raw matrix"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete number of dimensions of the raw matrix"""
            raise NotImplementedError("method is abstract")
        return locals()
    mat_ndim = property(**mat_ndim())

    def mat_shape():
        doc = "Shape of the raw matrix property."
        def fget(self):
            """Get the shape of the raw matrix"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the shape of the raw matrix"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the shape of the raw matrix"""
            raise NotImplementedError("method is abstract")
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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(self, values, axis, **kwargs):
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

    def delete(self, obj, axis, **kwargs):
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
        raise NotImplementedError("static method is abstract")

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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    def select(self, indices, axis, **kwargs):
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
        raise NotImplementedError("method is abstract")

    @staticmethod
    def concat(mats, axis, **kwargs):
        """
        Concatenate matrices together along an axis.

        Parameters
        ----------
        mats : array_like of matrices
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



################################################################################
################################## Utilities ###################################
################################################################################
def is_Matrix(v):
    """
    Determine whether an object is a Matrix.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a Matrix object instance.
    """
    return isinstance(v, Matrix)

def check_is_Matrix(v, varname):
    """
    Check if object is of type Matrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_Matrix(v):
        raise TypeError("'{0}' must be of type Matrix.".format(varname))

def cond_check_is_Matrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type Matrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a Matrix.
    """
    if cond(v):
        check_is_Matrix(v, varname)
