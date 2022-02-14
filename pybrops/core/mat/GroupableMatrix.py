from pybrops.core.mat.SortableMatrix import SortableMatrix

class GroupableMatrix(SortableMatrix):
    """
    An abstract class for groupable matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix in-place matrix grouping routines.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        GroupableMatrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(GroupableMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Grouping Methods ###################
    def group(self, axis, **kwargs):
        """
        Sort the Matrix, then populate grouping indices.

        Parameters
        ----------
        axis : int
            The axis along which values are grouped.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def is_grouped(self, axis, **kwargs):
        """
        Determine whether the Matrix has been sorted and grouped.

        Parameters
        ----------
        axis: int
            Axis along which to determine whether elements have been sorted and
            grouped.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        grouped : bool
            True or False indicating whether the Matrix has been sorted and
            grouped.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GroupableMatrix(v):
    """
    Determine whether an object is a GroupableMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a GroupableMatrix object instance.
    """
    return isinstance(v, GroupableMatrix)

def check_is_GroupableMatrix(v, varname):
    """
    Check if object is of type GroupableMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_GroupableMatrix(v):
        raise TypeError("'{0}' must be a GroupableMatrix".format(varname))

def cond_check_is_GroupableMatrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type GroupableMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        GroupableMatrix.
    """
    if cond(v):
        check_is_GroupableMatrix(v, varname)
