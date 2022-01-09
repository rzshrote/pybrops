class MatingProtocol:
    """docstring for MatingProtocol."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for abstract class MatingProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(MatingProtocol, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgmat, sel, ncross, nprogeny, miscout, **kwargs):
        """
        Mate individuals according to a mating scheme.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            A PhasedGenotypeMatrix of parental candidates.
        sel : numpy.ndarray
            Array of indices specifying a cross pattern. Each index corresponds
            to an individual in 'pgvmat'.
        ncross : numpy.ndarray
            Number of crosses to perform per cross pattern.
        nprogeny : numpy.ndarray
            Number of progeny to generate per cross.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : PhasedGenotypeMatrix
            A PhasedGenotypeMatrix of progeny.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_MatingProtocol(v):
    """
    Determine whether an object is a MatingProtocol.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a MatingProtocol object instance.
    """
    return isinstance(v, MatingProtocol)

def check_is_MatingProtocol(v, varname):
    """
    Check if object is of type MatingProtocol. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, MatingProtocol):
        raise TypeError("'%s' must be a MatingProtocol." % varname)

def cond_check_is_MatingProtocol(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type MatingProtocol. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a MatingProtocol.
    """
    if cond(v):
        check_is_MatingProtocol(v, varname)
