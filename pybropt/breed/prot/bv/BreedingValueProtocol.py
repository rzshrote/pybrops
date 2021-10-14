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
        **kwargs : dict
            Additional keyword arguments.
        """
        super(BreedingValueProtocol, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def estimate(self, obj, gmat, **kwargs):
        """
        Estimate breeding values.

        Parameters
        ----------
        obj : PhenotypeDataFrame, BreedingValueMatrix
            Phenotype dataframe or breeding value matrix to use to estimate
            breeding values.
        gmat : GenotypeMatrix
            Genotype matrix to use for estimation. Also used to align genotypes
            in estimation output.

        Returns
        -------
        bvmat : BreedingValueMatrix
            Breeding value matrix.
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
