from . import BreedingEdge

class EmigrationOperator(BreedingEdge):
    """docstring for EmigrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class EmigrationOperator.

        Parameters
        ----------
        **kwargs : **dict
            Additional keyword arguments.
        """
        super(EmigrationOperator, self).__init__(
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def emigrate(self, bnode, **kwargs):
        """
        Emigrate individuals to a BreedingNode.

        Parameters
        ----------
        bnode : BreedingNode
            A BreedingNode object to add individuals.
        kwargs : dict
            Dictionary of data for individuals to be added.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_EmigrationOperator(v):
    """
    Determine whether an object is a EmigrationOperator.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a EmigrationOperator object instance.
    """
    return isinstance(v, EmigrationOperator)

def check_is_EmigrationOperator(v, varname):
    """
    Check if object is of type EmigrationOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, EmigrationOperator):
        raise TypeError("'%s' must be a EmigrationOperator." % varname)

def cond_check_is_EmigrationOperator(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type EmigrationOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a EmigrationOperator.
    """
    if cond(v):
        check_is_EmigrationOperator(v, varname)
