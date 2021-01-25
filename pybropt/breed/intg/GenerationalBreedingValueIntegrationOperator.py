from . import BreedingValueIntegrationOperator

class GenerationalBreedingValueIntegrationOperator(BreedingValueIntegrationOperator):
    """docstring for GenerationalBreedingValueIntegrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, gqlen, gwind, gmult, **kwargs):
        super(GenerationalBreedingValueIntegrationOperator, self).__init__(**kwargs)
        self.gqlen = gqlen
        self.gwind = gwind
        self.gmult = gmult

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

        # calculate the taxa_grp minimum threshold
        taxa_min = (t_cur - (self.gqlen + self.gwind)) * self.gmult

        # process breeding value matrix
        bval_new["main"] = bvmat
        # mask = bval_new["main"].taxa_grp < taxa_min     # create breeding value mask
        # bval_new["main"].delete(mask, axis = 0)         # delete old taxa
        # bval_new["main"].append(                        # add new taxa
        #     values = bvmat.mat,
        #     axis = 0,
        #     raw = bvmat.raw,    # FIXME: will raise attribute error if not a DenseEstimatedBreedingValueMatrix
        #     taxa = bvmat.taxa,
        #     taxa_grp = bvmat.taxa_grp
        # )

        # process true breeding value matrix
        bval_new["main_true"] = bvmat_true
        # mask = bval_new["main_true"].taxa_grp < taxa_min     # create breeding value mask
        # bval_new["main_true"].delete(mask, axis = 0)         # delete old taxa
        # bval_new["main_true"].append(                        # add new taxa
        #     values = bvmat_true.mat,
        #     axis = 0,
        #     taxa = bvmat_true.taxa,
        #     taxa_grp = bvmat_true.taxa_grp
        # )

        misc = {}

        return bval_new, misc



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenerationalBreedingValueIntegrationOperator(v):
    return isinstance(v, GenerationalBreedingValueIntegrationOperator)

def check_is_GenerationalBreedingValueIntegrationOperator(v, vname):
    if not isinstance(v, GenerationalBreedingValueIntegrationOperator):
        raise TypeError("variable '{0}' must be a GenerationalBreedingValueIntegrationOperator".format(vname))

def cond_check_is_GenerationalBreedingValueIntegrationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GenerationalBreedingValueIntegrationOperator(v, vname)
