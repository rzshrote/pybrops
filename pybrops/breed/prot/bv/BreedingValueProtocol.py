class BreedingValueProtocol:
    """docstring for BreedingValueProtocol."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class BreedingValueProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(BreedingValueProtocol, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def estimate(self, ptobj, gtobj, miscout, **kwargs):
        """
        Estimate breeding values.

        Parameters
        ----------
        ptobj : PhenotypeDataFrame
            An object containing phenotype data. Must be a phenotype data frame.
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a genotype matrix.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : BreedingValueMatrix
            A matrix of breeding values.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingValueProtocol(v):
    """
    Determine whether an object is a BreedingValueProtocol.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a BreedingValueProtocol object instance.
    """
    return isinstance(v, BreedingValueProtocol)

def check_is_BreedingValueProtocol(v, varname):
    """
    Check if object is of type BreedingValueProtocol. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingValueProtocol):
        raise TypeError("'%s' must be a BreedingValueProtocol." % varname)

def cond_check_is_BreedingValueProtocol(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type BreedingValueProtocol. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a BreedingValueProtocol.
    """
    if cond(v):
        check_is_BreedingValueProtocol(v, varname)