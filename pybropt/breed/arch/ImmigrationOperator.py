from pybropt.breed.arch.BreedingEdge import BreedingEdge

class ImmigrationOperator(BreedingEdge):
    """docstring for ImmigrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class ImmigrationOperator.

        Parameters
        ----------
        **kwargs : **dict
            Additional keyword arguments.
        """
        super(ImmigrationOperator, self).__init__(
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def immigrate(self, bnode, **kwargs):
        """
        Immigrate individuals from a BreedingNode.

        Parameters
        ----------
        bnode : BreedingNode
            A BreedingNode object from which to pull individuals.

        Returns
        -------
        kwindiv : dict
            Dictionary of data from individuals selected.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_ImmigrationOperator(v):
    """
    Determine whether an object is a ImmigrationOperator.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a ImmigrationOperator object instance.
    """
    return isinstance(v, ImmigrationOperator)

def check_is_ImmigrationOperator(v, varname):
    """
    Check if object is of type ImmigrationOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, ImmigrationOperator):
        raise TypeError("'%s' must be a ImmigrationOperator." % varname)

def cond_check_is_ImmigrationOperator(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type ImmigrationOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a ImmigrationOperator.
    """
    if cond(v):
        check_is_ImmigrationOperator(v, varname)
