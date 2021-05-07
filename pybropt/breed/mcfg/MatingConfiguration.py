class MatingConfiguration:
    """docstring for MatingConfiguration."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(MatingConfiguration, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def config(self):
        """
        Construct a mating configuration from information within self.

        Returns
        -------
        cfg : dict
            A dictionary specifying the mating configuration.
            Dictionary must have the following fields:
                pgvmat : PhasedGenotypeVariantMatrix
                    A GenotypeVariantMatrix containing candidate breeding individuals.
                sel : numpy.ndarray
                    An array of indices indicating the cross pattern.
                ncross : numpy.ndarray
                    Number of cross patterns to perform.
                nprogeny : numpy.ndarray
                    Number of progeny to generate per cross.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_MatingConfiguration(v):
    return isinstance(v, MatingConfiguration)

def check_is_MatingConfiguration(v, varname):
    if not isinstance(v, MatingConfiguration):
        raise TypeError("'%s' must be a MatingConfiguration." % varname)

def cond_check_is_MatingConfiguration(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_MatingConfiguration(v, varname)
