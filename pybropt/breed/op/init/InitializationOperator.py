class InitializationOperator:
    """docstring for InitializationOperator."""

    def __init__(self, **kwargs):
        """
        Constructor for the abstract class InitializationOperator.

        Parameters
        ----------
        **kwargs : **dict
            Additional keyword arguments.
        """
        super(InitializationOperator, self).__init__()

    def initialize(self, miscout, **kwargs):
        """
        Initialize a breeding program.

        Parameters
        ----------
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple of length 5: (genome, geno, pheno, bval, gmod)
            Where:
                genome : dict
                    A dictionary of genomes for the breeding program.
                geno : dict
                    A dictionary of genotypes for the breeding program.
                pheno : dict
                    A dictionary of phenotypes for the breeding program.
                bval : dict
                    A dictionary of breeding values for the breeding program.
                gmod : dict
                    A dictionary of genomic models for the breeding program.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_InitializationOperator(v):
    """
    Determine whether an object is a InitializationOperator.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a InitializationOperator object instance.
    """
    return isinstance(v, InitializationOperator)

def check_is_InitializationOperator(v, varname):
    """
    Check if object is of type InitializationOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, InitializationOperator):
        raise TypeError("'%s' must be a InitializationOperator." % varname)

def cond_check_is_InitializationOperator(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type InitializationOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a InitializationOperator.
    """
    if cond(v):
        check_is_InitializationOperator(v, varname)
