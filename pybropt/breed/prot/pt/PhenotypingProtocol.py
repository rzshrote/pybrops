class PhenotypingProtocol:
    """docstring for PhenotypingProtocol."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(PhenotypingProtocol, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############### Genomic Model Properties ###############
    def gpmod():
        doc = "Genomic prediction model."
        def fget(self):
            """Get genomic prediction model"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set genomic prediction model"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete genomic prediction model"""
            raise NotImplementedError("method is abstract")
        return locals()
    gpmod = property(**gpmod())

    ################ Stochastic Parameters #################
    def var_E():
        doc = "Environmental variance for each trait."
        def fget(self):
            """Get environmental variance"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set environmental variance"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete environmental variance"""
            raise NotImplementedError("method is abstract")
        return locals()
    var_E = property(**var_E())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def phenotype(self, pgmat, **kwargs):
        """
        Phenotype a set of genotypes using a genomic prediction model.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes of the individuals to phenotype.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhenotypeDataFrame
            DataFrame containing phenotypes.
        """
        raise NotImplementedError("method is abstract")

    def set_h2(self, h2, pgmat, **kwargs):
        """
        Set the narrow sense heritability for environments.

        Parameters
        ----------
        h2 : float, numpy.ndarray
            Narrow sense heritability.
        pgmat : PhasedGenotypeVariantMatrix
            Founder genotypes.
        **kwargs : dict
            Additional keyword arguments
        """
        raise NotImplementedError("method is abstract")

    def set_H2(self, H2, pgmat, **kwargs):
        """
        Set the broad sense heritability for environments.

        Parameters
        ----------
        H2 : float, numpy.ndarray
            Broad sense heritability.
        pgmat : PhasedGenotypeVariantMatrix
            Founder genotypes.
        **kwargs : dict
            Additional keyword arguments
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhenotypingProtocol(v):
    return isinstance(v, PhenotypingProtocol)

def check_is_PhenotypingProtocol(v, vname):
    if not isinstance(v, PhenotypingProtocol):
        raise TypeError("variable '{0}' must be a PhenotypingProtocol".format(vname))

def cond_check_is_PhenotypingProtocol(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_PhenotypingProtocol(v, vname)
