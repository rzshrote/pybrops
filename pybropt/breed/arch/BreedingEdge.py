class BreedingEdge:
    """docstring for BreedingEdge."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class BreedingEdge.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(BreedingEdge, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingEdge(v):
    """
    Determine whether an object is a BreedingEdge.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a BreedingEdge object instance.
    """
    return isinstance(v, BreedingEdge)

def check_is_BreedingEdge(v, varname):
    """
    Check if object is of type BreedingEdge. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingEdge):
        raise TypeError("'%s' must be a BreedingEdge." % varname)

def cond_check_is_BreedingEdge(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type BreedingEdge. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a BreedingEdge.
    """
    if cond(v):
        check_is_BreedingEdge(v, varname)
