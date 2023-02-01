from . import BreedingValueIntegrationOperator

class SimpleBreedingValueIntegrationOperator(BreedingValueIntegrationOperator):
    """docstring for SimpleBreedingValueIntegrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict):
        super(SimpleBreedingValueIntegrationOperator, self).__init__(**kwargs)

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
        # copy dictionaries
        bval_new = dict(bval)

        # process breeding value matrix
        bval_new["main"] = bvmat

        # process true breeding value matrix
        bval_new["main_true"] = bvmat_true

        # empty dict
        misc = {}

        return bval_new, misc



################################################################################
################################## Utilities ###################################
################################################################################
def is_SimpleBreedingValueIntegrationOperator(v):
    return isinstance(v, SimpleBreedingValueIntegrationOperator)

def check_is_SimpleBreedingValueIntegrationOperator(v, vname):
    if not isinstance(v, SimpleBreedingValueIntegrationOperator):
        raise TypeError("variable '{0}' must be a SimpleBreedingValueIntegrationOperator".format(vname))

def cond_check_is_SimpleBreedingValueIntegrationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_SimpleBreedingValueIntegrationOperator(v, vname)
