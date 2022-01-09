class BreedingGraph:
    """docstring for BreedingGraph."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class BreedingGraph.

        Parameters
        ----------
        **kwargs : **dict
            Additional keyword arguments.
        """
        super(BreedingGraph, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################### Graph properties ###################
    def graph():
        doc = "Access to the underlying graph structure."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    graph = property(**graph())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingGraph(v):
    """
    Determine whether an object is a BreedingGraph.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a BreedingGraph object instance.
    """
    return isinstance(v, BreedingGraph)

def check_is_BreedingGraph(v, varname):
    """
    Check if object is of type BreedingGraph. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingGraph):
        raise TypeError("'%s' must be a BreedingGraph." % varname)

def cond_check_is_BreedingGraph(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type BreedingGraph. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a BreedingGraph.
    """
    if cond(v):
        check_is_BreedingGraph(v, varname)
