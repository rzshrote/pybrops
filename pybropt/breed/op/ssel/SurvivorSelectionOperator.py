class SurvivorSelectionOperator:
    """docstring for SurvivorSelectionOperator."""

    def __init__(self, **kwargs):
        """
        Constructor for the abstract class SurvivorSelectionOperator.

        Parameters
        ----------
        **kwargs : **dict
            Additional keyword arguments.
        """
        super(SurvivorSelectionOperator, self).__init__()

    def sselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        """
        Select progeny survivors in a breeding program.

        Parameters
        ----------
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
def is_SurvivorSelectionOperator(v):
    """
    Determine whether an object is a SurvivorSelectionOperator.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a SurvivorSelectionOperator object instance.
    """
    return isinstance(v, SurvivorSelectionOperator)

def check_is_SurvivorSelectionOperator(v, varname):
    """
    Check if object is of type SurvivorSelectionOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SurvivorSelectionOperator):
        raise TypeError("'%s' must be a SurvivorSelectionOperator." % varname)

def cond_check_is_SurvivorSelectionOperator(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type SurvivorSelectionOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a SurvivorSelectionOperator.
    """
    if cond(v):
        check_is_SurvivorSelectionOperator(v, varname)
