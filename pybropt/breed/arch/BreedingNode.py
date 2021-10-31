class BreedingNode:
    """docstring for BreedingNode."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class BreedingNode.

        Parameters
        ----------
        **kwargs : **dict
            Additional keyword arguments.
        """
        super(BreedingNode, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############# Generation number properties #############
    def t_cur():
        doc = "Current generation number of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    t_cur = property(**t_cur())

    def t_max():
        doc = "Maximum generation number of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    t_max = property(**t_max())

    ################ Population properties #################
    def geno():
        doc = "Main breeding population of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    geno = property(**geno())

    ############## Breeding value properties ###############
    def bval():
        doc = "Dictionary of breeding values for the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    bval = property(**bval())

    ############### Genomic model properties ###############
    def gmod():
        doc = "Dictionary of genomic models for the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    gmod = property(**gmod())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingNode(v):
    """
    Determine whether an object is a BreedingNode.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a BreedingNode object instance.
    """
    return isinstance(v, BreedingNode)

def check_is_BreedingNode(v, varname):
    """
    Check if object is of type BreedingNode. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingNode):
        raise TypeError("'%s' must be a BreedingNode." % varname)

def cond_check_is_BreedingNode(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type BreedingNode. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a BreedingNode.
    """
    if cond(v):
        check_is_BreedingNode(v, varname)
