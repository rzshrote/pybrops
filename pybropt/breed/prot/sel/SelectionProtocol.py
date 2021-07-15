class SelectionProtocol:
    """docstring for SelectionProtocol."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(SelectionProtocol, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
        """
        Select parents individuals for breeding.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing five objects: (pgmat, sel, ncross, nprogeny, misc)
            pgmat : PhasedGenotypeMatrix
                A PhasedGenotypeMatrix of parental candidates.
            sel : numpy.ndarray
                Array of indices specifying a cross pattern. Each index
                corresponds to an individual in 'pgmat'.
            ncross : numpy.ndarray
                Number of crosses to perform per cross pattern.
            nprogeny : numpy.ndarray
                Number of progeny to generate per cross.
            misc : dict
                Miscellaneous output (user defined).
        """
        raise NotImplementedError("method is abstract")

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
        """
        Return an objective function constructed from inputs for fast calling.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : function
            Output function
        """
        raise NotImplementedError("method is abstract")

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
        """
        Return a vectorized objective function constructed from inputs for fast
        calling.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : function
            Output function
        """
        raise NotImplementedError("method is abstract")

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
        """
        Calculate a Pareto frontier for objectives.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing three objects (frontier, sel_config, misc)
            Elements
            --------
            frontier : numpy.ndarray
                Array of shape (q,v) containing Pareto frontier points.
                Where:
                    'q' is the number of points in the frontier.
                    'v' is the number of objectives for the frontier.
            sel_config : numpy.ndarray
                Array of shape (q,k) containing parent selection decisions for
                each corresponding point in the Pareto frontier.
                Where:
                    'q' is the number of points in the frontier.
                    'k' is the number of search space decision variables.
            misc : dict
                A dictionary of miscellaneous output. (User specified)
        """
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, *args):
        """
        Objective function for the selection protocol.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A parent contribution vector of shape (n,) and floating dtype.
            Where:
                'n' is the number of individuals.
        *args : tuple
            Additional arguments.

        Returns
        -------
        out : numpy.ndarray, scalar
            An array of objective function values of shape (o,) or a scalar.
            Where:
                'o' is the number of objectives.
        """
        raise NotImplementedError("method is abstract")

    @staticmethod
    def objfn_vec_static(sel, *args):
        """
        Vectorized objective function for the selection protocol.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A parent contribution vector of shape (j,n) and floating dtype.
            Where:
                'j' is the number of selection configurations.
                'n' is the number of individuals.
        *args : tuple
            Additional arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of objective function values of shape (j,o) or (j,).
            Where:
                'j' is the number of selection configurations.
                'o' is the number of objectives.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_SelectionProtocol(v):
    return isinstance(v, SelectionProtocol)

def check_is_SelectionProtocol(v, vname):
    if not isinstance(v, SelectionProtocol):
        raise TypeError("variable '{0}' must be a SelectionProtocol".format(vname))

def cond_check_is_SelectionProtocol(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_SelectionProtocol(v, vname)
