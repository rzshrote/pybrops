class MatingProtocol:
    """docstring for MatingProtocol."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(MatingProtocol, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def mate(self, pgmat, sel, ncross, nprogeny, **kwargs):
        """
        Mate individuals according to a mating scheme.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            A PhasedGenotypeMatrix of parental candidates.
        sel : numpy.ndarray
            Array of indices specifying a cross pattern. Each index corresponds
            to an individual in 'pgvmat'.
        ncross : numpy.ndarray
            Number of crosses to perform per cross pattern.
        nprogeny : numpy.ndarray
            Number of progeny to generate per cross.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing two elements: (progeny, misc)
            progeny : PhasedGenotypeMatrix
                A PhasedGenotypeMatrix of progeny.
            misc : dict
                Miscellaneous output (user defined).
        """
        raise NotImplementedError("method is abstract")
