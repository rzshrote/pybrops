class TrueEvaluationOperator(EvaluationOperator):
    """docstring for TrueEvaluationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(TrueEvaluationOperator, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def evaluate(self, t_cur, t_max, pgvmat, gmod, **kwargs):
        """
        Parameters
        ----------
        pgvmat : PhasedGenotypeVariantMatrix
        gmod : dict

        Returns
        -------
        out : tuple
            A tuple containing two elements: (bvmat, misc)
            bvmat : BreedingValueMatrix
                A matrix of breeding values
            bvmat_true : BreedingValueMatrix
                A matrix of true breeding values
            misc : dict
                Miscellaneous output (user defined).
        """
        bvmat = gmod["true"].pred(pgvmat)   # make predictions using true model
        misc = {}                           # declare an empty dict

        return bvmat, bvmat, misc



################################################################################
################################## Utilities ###################################
################################################################################
def is_TrueEvaluationOperator(v):
    return isinstance(v, TrueEvaluationOperator)

def check_is_TrueEvaluationOperator(v, vname):
    if not isinstance(v, TrueEvaluationOperator):
        raise TypeError("variable '{0}' must be a TrueEvaluationOperator".format(vname))

def cond_check_is_TrueEvaluationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_TrueEvaluationOperator(v, vname)
