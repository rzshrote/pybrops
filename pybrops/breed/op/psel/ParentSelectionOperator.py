"""
Module defining interfaces and associated error checking routines for
breeding program parental selection operators.
"""

class ParentSelectionOperator:
    """
    Abstract class defining interfaces for parental selection within an entire
    breeding program.

    The purpose of this abstract class is to provide functionality for:
        1) Parental selection for an entire breeding program.
    """

    def __init__(self, **kwargs):
        """
        Constructor for the abstract class ParentSelectionOperator.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(ParentSelectionOperator, self).__init__()

    def pselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, miscout, **kwargs):
        """
        Select individuals to serve as parents in a breeding program.

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
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple of length 6: ``(mcfg, genome, geno, pheno, bval, gmod)``

            Where:

            - ``mcfg`` is a ``dict`` of mating configurations for the breeding program.
            - ``genome`` is a ``dict`` of genomes for the breeding program.
            - ``geno`` is a ``dict`` of genotypes for the breeding program.
            - ``pheno`` is a ``dict`` of phenotypes for the breeding program.
            - ``bval`` is a ``dict`` of breeding values for the breeding program.
            - ``gmod`` is a ``dict`` of genomic models for the breeding program.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_ParentSelectionOperator(v):
    """
    Determine whether an object is a ParentSelectionOperator.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a ParentSelectionOperator object instance.
    """
    return isinstance(v, ParentSelectionOperator)

def check_is_ParentSelectionOperator(v, varname):
    """
    Check if object is of type ParentSelectionOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, ParentSelectionOperator):
        raise TypeError("'%s' must be a ParentSelectionOperator." % varname)

def cond_check_is_ParentSelectionOperator(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type ParentSelectionOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a ParentSelectionOperator.
    """
    if cond(v):
        check_is_ParentSelectionOperator(v, varname)
