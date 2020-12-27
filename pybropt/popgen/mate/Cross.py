class Cross:
    """docstring for Cross."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(Cross, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgvmat, sel, ncross, nprogeny, **kwargs):
        """
        Mate individuals according to a mate selection scheme.

        Parameters
        ----------
        pgvmat : PhasedGenotypeVariantMatrix
            A GenotypeVariantMatrix from which to mate.
        sel : numpy.ndarray
            Cross pattern.
        ncross : int, numpy.ndarray
            Number of cross patterns to perform.
        nprogeny : int, numpy.ndarray
            Number of progeny to generate per cross.

        Returns
        -------
        progeny : PhasedGenotypeVariantMatrix
            A PhasedGenotypeVariantMatrix of progeny.
        """
        raise NotImplementedError("method is abstract")

    def mate_cfg(self, mcfg):
        """
        Mate individuals according to a MatingConfiguration.

        Parameters
        ----------
        mcfg : MatingConfiguration
            A MatingConfiguration from which to implement matings.

        Returns
        -------
        progeny : PhasedGenotypeVariantMatrix
            A PhasedGenotypeVariantMatrix of progeny.
        """
        raise NotImplementedError("method is abstract")


################################################################################
################################## Utilities ###################################
################################################################################
def is_Cross(v):
    return isinstance(v, Cross)

def check_is_Cross(v, varname):
    if not isinstance(v, Cross):
        raise TypeError("'%s' must be a Cross." % varname)

def cond_check_is_Cross(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_Cross(v, varname)
