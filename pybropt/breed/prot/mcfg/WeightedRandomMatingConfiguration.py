from . import MatingConfiguration

class WeightedRandomMatingConfiguration(MatingConfiguration):
    """docstring for WeightedRandomMatingConfiguration."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, arg):
        super(WeightedRandomMatingConfiguration, self).__init__()
        self.arg = arg

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
                pgvmat : PhasedGenotypeMatrix
                    A GenotypeMatrix containing candidate breeding individuals.
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
def is_WeightedRandomMatingConfiguration(v):
    return isinstance(v, WeightedRandomMatingConfiguration)

def check_is_WeightedRandomMatingConfiguration(v, varname):
    if not isinstance(v, WeightedRandomMatingConfiguration):
        raise TypeError("'%s' must be a WeightedRandomMatingConfiguration." % varname)

def cond_check_is_WeightedRandomMatingConfiguration(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_WeightedRandomMatingConfiguration(v, varname)
