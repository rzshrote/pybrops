class MatingOperator:
    """docstring for MatingOperator."""

    def __init__(self, **kwargs):
        """
        Constructor for the abstract class MatingOperator

        Parameters
        ----------
        **kwargs : **dict
            Additional keyword arguments.
        """
        super(MatingOperator, self).__init__()

    def mate(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        """
        Mate individuals selected as parents in a breeding program.

        Parameters
        ----------
        mcfg : dict
            Dictionary of mating configurations for the breeding program.
        genome : dict
            Dictionary of genomes for the breeding program.
        geno : dict
            Dictionary of genotypes for the breeding program.
        pheno : dict
            Dictionary of phenotypes for the breeding program.
        bval : dict
            Dictionary of breeding values for the breeding program.
        gmod : dict
            Dictionary of genomic models for the breeding program.
        t_cur : int
            Current time in the breeding program.
        t_max : int
            Deadline time for the breeding program.
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
def is_MatingOperator(v):
    """
    Determine whether an object is a MatingOperator.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a MatingOperator object instance.
    """
    return isinstance(v, MatingOperator)

def check_is_MatingOperator(v, varname):
    """
    Check if object is of type MatingOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, MatingOperator):
        raise TypeError("'%s' must be a MatingOperator." % varname)

def cond_check_is_MatingOperator(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type MatingOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a MatingOperator.
    """
    if cond(v):
        check_is_MatingOperator(v, varname)
