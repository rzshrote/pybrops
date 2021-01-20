class BreedingValueIntegrationOperator:
    """docstring for BreedingValueIntegrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(BreedingValueIntegrationOperator, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def bvintegrate(self, t_cur, t_max, bvmat, bvmat_true, bval):
        """
        Integrate breeding values into bval dictionary.

        Parameters
        ----------
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        bvmat : BreedingValueMatrix
            Estimated breeding values
        bvmat_true : BreedingValueMatrix
            True breeding values
        bval : dict
            Breeding value dictionary into which to integrate.

        Returns
        -------
        out : tuple
            (bval_new, misc)
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingValueIntegrationOperator(v):
    return isinstance(v, BreedingValueIntegrationOperator)

def check_is_BreedingValueIntegrationOperator(v, vname):
    if not isinstance(v, BreedingValueIntegrationOperator):
        raise TypeError("variable '{0}' must be a BreedingValueIntegrationOperator".format(vname))

def cond_check_is_BreedingValueIntegrationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_BreedingValueIntegrationOperator(v, vname)
