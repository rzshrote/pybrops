from pybropt.breed.arch.BreedingNode import BreedingNode

class BreedingProgram(BreedingNode):
    """docstring for BreedingProgram."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class BreedingProgram.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(BreedingProgram, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############ Starting condition containers #############
    def start_genome():
        doc = "Starting genomes for individuals in the breeding program."
        def fget(self):
            """Get starting genomes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set starting genomes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete starting genomes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    start_genome = property(**start_genome())

    def start_geno():
        doc = "Starting genotypes for individuals in the breeding program."
        def fget(self):
            """Get starting genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set starting genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete starting genotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    start_geno = property(**start_geno())

    def start_pheno():
        doc = "Starting phenotypes for individuals in the breeding program."
        def fget(self):
            """Get starting phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set starting phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete starting phenotypes for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    start_pheno = property(**start_pheno())

    def start_bval():
        doc = "Starting breeding values for individuals in the breeding program."
        def fget(self):
            """Get starting breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set starting breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete starting breeding values for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    start_bval = property(**start_bval())

    def start_gmod():
        doc = "Starting genomic models for individuals in the breeding program."
        def fget(self):
            """Get starting genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set starting genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete starting genomic models for individuals in the breeding program"""
            raise NotImplementedError("method is abstract")
        return locals()
    start_gmod = property(**start_gmod())

    ######### Breeding program operator properties #########
    def initop():
        doc = "Initialization operator"
        def fget(self):
            """Get the initialization operator"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the initialization operator"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the initialization operator"""
            raise NotImplementedError("method is abstract")
        return locals()
    initop = property(**initop())

    def pselop():
        doc = "Parent selection operator."
        def fget(self):
            """Get the parent selection operator"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the parent selection operator"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the parent selection operator"""
            raise NotImplementedError("method is abstract")
        return locals()
    pselop = property(**pselop())

    def mateop():
        doc = "Mating operator."
        def fget(self):
            """Get the mating operator"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the mating operator"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the mating operator"""
            raise NotImplementedError("method is abstract")
        return locals()
    mateop = property(**mateop())

    def evalop():
        doc = "Evaluation operator."
        def fget(self):
            """Get the evaluation operator"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the evaluation operator"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the evaluation operator"""
            raise NotImplementedError("method is abstract")
        return locals()
    evalop = property(**evalop())

    def sselop():
        doc = "Survivor selection operator."
        def fget(self):
            """Get the survivor selection operator"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the survivor selection operator"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the survivor selection operator"""
            raise NotImplementedError("method is abstract")
        return locals()
    sselop = property(**sselop())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############# Initialize breeding program ##############
    def initialize(self, **kwargs):
        """
        Initialize the breeding program with genotypes, phenotypes, and genomic
        models.
        """
        raise NotImplementedError("method is abstract")

    def is_initialized(self, **kwargs):
        """
        Return whether or not the BreedingProgram has been initialized with a
        starting set of conditions.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : boolean
            True if the BreedingProgram has been initialized.
            False if the BreedingProgram has not been initialized.
        """
        raise NotImplementedError("method is abstract")

    ############# Population evolution methods #############
    def reset(self, **kwargs):
        """
        Reset the evolution of the breeding program back to starting conditions.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def advance(self, ngen, lbook, **kwargs):
        """
        Advance the breeding program by a specified number of generations.

        Parameters
        ----------
        ngen : int
            Number of generations to advance the BreedingProgram.
        lbook : Logbook
            Logbook into which to write statistics.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def evolve(self, nrep, ngen, lbook, **kwargs):
        """
        Evolve the breeding program for a set number of replications and
        generations. The BreedingProgram is restarted using the starting geno,
        bval, gmod containers.

        Parameters
        ----------
        nrep : int
            Number of evolution replicates.
        ngen : int
            Number of generations to evolve the population for each replicate.
            Note that this does not modify 't_max'.
        lbook : Logbook
            Logbook into which to write statistics.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingProgram(v):
    """
    Determine whether an object is a BreedingProgram.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a BreedingProgram object instance.
    """
    return isinstance(v, BreedingProgram)

def check_is_BreedingProgram(v, varname):
    """
    Check if object is of type BreedingProgram. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingProgram):
        raise TypeError("'%s' must be a BreedingProgram." % varname)

def cond_check_is_BreedingProgram(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type BreedingProgram. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a BreedingProgram.
    """
    if cond(v):
        check_is_BreedingProgram(v, varname)
