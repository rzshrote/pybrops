from pybrops.core.mat.SortableMatrix import SortableMatrix

class TraitMatrix(SortableMatrix):
    """
    An abstract class for matrix wrapper objects with trait metadata.

    The purpose of this abstract class is to provide base functionality for:
        1) trait metadata manipulation routines.
        2) trait manipulation routines.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        TraitMatrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(TraitMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ###################### Trait data ######################
    def trait():
        doc = "Trait label property."
        def fget(self):
            """Get trait label array"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set trait label array"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete trait label array"""
            raise NotImplementedError("method is abstract")
        return locals()
    trait = property(**trait())

    #################### Trait metadata ####################
    def ntrait():
        doc = "Number of traits property."
        def fget(self):
            """Get number of traits"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set number of traits"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete number of traits"""
            raise NotImplementedError("method is abstract")
        return locals()
    ntrait = property(**ntrait())

    def trait_axis():
        doc = "Axis along which traits are stored property."
        def fget(self):
            """Get trait axis number"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set trait axis number"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete trait axis number"""
            raise NotImplementedError("method is abstract")
        return locals()
    trait_axis = property(**trait_axis())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin_trait(self, values, trait, **kwargs):
        """
        Add additional elements to the end of the Matrix along the trait axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        trait : numpy.ndarray
            Taxa names to adjoin to the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    def delete_trait(self, obj, **kwargs):
        """
        Delete sub-arrays along the trait axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    def insert_trait(self, obj, values, trait, **kwargs):
        """
        Insert values along the trait axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        trait : numpy.ndarray
            Taxa names to insert into the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    def select_trait(self, indices, **kwargs):
        """
        Select certain values from the Matrix along the trait axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            The output Matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        raise NotImplementedError("method is abstract")

    @staticmethod
    def concat_trait(mats, **kwargs):
        """
        Concatenate list of Matrix together along the trait axis.

        Parameters
        ----------
        mats : array_like of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        raise NotImplementedError("static method is abstract")

    ######### Matrix element in-place-manipulation #########
    def append_trait(self, values, trait, **kwargs):
        """
        Append values to the Matrix along the trait axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        trait : numpy.ndarray
            Taxa names to append to the Matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def remove_trait(self, obj, **kwargs):
        """
        Remove sub-arrays along the trait axis.

        Parameters
        ----------
        obj : slice, int, or array of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def incorp_trait(self, obj, values, trait, **kwargs):
        """
        Incorporate values along the trait axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        trait : numpy.ndarray
            Taxa names to incorporate into the Matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ################### Sorting Methods ####################
    def lexsort_trait(self, keys, **kwargs):
        """
        Perform an indirect stable sort using a sequence of keys along the trait
        axis.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        indices : (N,) ndarray of ints
            Array of indices that sort the keys along the specified axis.
        """
        raise NotImplementedError("method is abstract")

    def reorder_trait(self, indices, **kwargs):
        """
        Reorder elements of the Matrix along the trait axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def sort_trait(self, keys, **kwargs):
        """
        Sort slements of the Matrix along the trait axis using a sequence of
        keys. Note this modifies the Matrix in-place.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_TraitMatrix(v):
    """
    Determine whether an object is a TraitMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a TraitMatrix object instance.
    """
    return isinstance(v, TraitMatrix)

def check_is_TraitMatrix(v, varname):
    """
    Check if object is of type TraitMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_TraitMatrix(v):
        raise TypeError("'{0}' must be a TraitMatrix".format(varname))

def cond_check_is_TraitMatrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type TraitMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        TraitMatrix.
    """
    if cond(v):
        check_is_TraitMatrix(v, varname)
