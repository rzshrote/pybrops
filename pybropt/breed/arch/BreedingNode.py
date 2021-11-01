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

    ############ Program information containers ############
    def genome():
        doc = "Genomes for individuals in the breeding program."
        def fget(self):
            """Get genomes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set genomes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete genomes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    genome = property(**genome())

    def geno():
        doc = "Genotypes for individuals in the breeding program."
        def fget(self):
            """Get genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    geno = property(**geno())

    def pheno():
        doc = "Phenotypes for individuals in the breeding program."
        def fget(self):
            """Get phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    pheno = property(**pheno())

    def bval():
        doc = "Breeding values for individuals in the breeding program."
        def fget(self):
            """Get breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    bval = property(**bval())

    def gmod():
        doc = "Genomic models for individuals in the breeding program."
        def fget(self):
            """Get genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    gmod = property(**gmod())

    ############# Generation number properties #############
    def t_cur():
        doc = "Current time of the BreedingNode."
        def fget(self):
            """Get the current time for the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the current time for the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the current time for the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    t_cur = property(**t_cur())

    def t_max():
        doc = "Maximum time of the BreedingNode."
        def fget(self):
            """Get the maximum time for the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the maximum time for the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the maximum time for the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    t_max = property(**t_max())

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
