class Logbook:
    """docstring for Logbook."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class Logbook.

        Parameters
        ----------
        **kwargs : **dict
            Additional keyword arguments.
        """
        super(Logbook, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def data():
        doc = "The data property."
        def fget(self):
            """Get Logbook data"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set Logbook data"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete Logbook data"""
            raise NotImplementedError("method is abstract")
        return locals()
    data = property(**data())

    def rep():
        doc = "The rep property."
        def fget(self):
            """Get replicate number"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set replicate number"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete replicate number"""
            raise NotImplementedError("method is abstract")
        return locals()
    rep = property(**rep())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def log_initialize(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        """
        Record information directly after 'InitializationOperator.initialize'
        is called.

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
        **kwargs : **dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def log_pselect(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        """
        Record information directly after 'ParentSelectionOperator.pselect'
        is called.

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
        **kwargs : **dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def log_mate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        """
        Record information directly after 'MatingOperator.mate' is called.

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
        **kwargs : **dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def log_evaluate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        """
        Record information directly after 'EvaluationOperator.evaluate' is
        called.

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
        **kwargs : **dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def log_sselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs):
        """
        Record information directly after 'SurvivorSelectionOperator.sselect'
        is called.

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
        **kwargs : **dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def reset(self):
        """
        Reset Logbook internals.
        """
        raise NotImplementedError("method is abstract")

    def write(self, fname):
        """
        Write Logbook to file

        Parameters
        ----------
        fname : str
            File name to which to write file.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_Logbook(v):
    """
    Determine whether an object is a Logbook.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a Logbook object instance.
    """
    return isinstance(v, Logbook)

def check_is_Logbook(v, varname):
    """
    Check if object is of type Logbook. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, Logbook):
        raise TypeError("'%s' must be a Logbook." % varname)

def cond_check_is_Logbook(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type Logbook. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a Logbook.
    """
    if cond(v):
        check_is_Logbook(v, varname)
