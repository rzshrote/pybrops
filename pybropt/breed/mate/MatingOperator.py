class MatingOperator:
    """docstring for MatingOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(MatingOperator, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, t_cur, t_max, pgvmat, sel, ncross, nprogeny, rng, **kwargs):
        """
        Mate individuals according to a mate selection scheme.

        Parameters
        ----------
        pgvmat : PhasedGenotypeVariantMatrix
            A PhasedGenotypeVariantMatrix of parental candidates.
        sel : numpy.ndarray
            Array of indices specifying a cross pattern. Each index corresponds
            to an individual in 'pgvmat'.
        ncross : numpy.ndarray
            Number of crosses to perform per cross pattern.
        nprogeny : numpy.ndarray
            Number of progeny to generate per cross.
        rng : random number generator instance
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing two elements: (progeny, misc)
            progeny : PhasedGenotypeVariantMatrix
                A PhasedGenotypeVariantMatrix of progeny.
            misc : dict
                Miscellaneous output (user defined).
        """
        raise NotImplementedError("method is abstract")


################################################################################
################################## Utilities ###################################
################################################################################
def is_MatingOperator(v):
    return isinstance(v, MatingOperator)

def check_is_MatingOperator(v, vname):
    if not isinstance(v, MatingOperator):
        raise TypeError("variable '{0}' must be a MatingOperator".format(vname))

def cond_check_is_MatingOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_MatingOperator(v, vname)
