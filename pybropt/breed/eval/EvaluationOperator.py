class EvaluationOperator:
    """docstring for EvaluationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(EvaluationOperator, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def evaluate(self, pgvmat, gmod, **kwargs):
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
            misc : dict
                Miscellaneous output (user defined).
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_EvaluationOperator(v):
    return isinstance(v, EvaluationOperator)

def check_is_EvaluationOperator(v, varname):
    if not isinstance(v, EvaluationOperator):
        raise TypeError("'%s' must be a EvaluationOperator." % varname)

def cond_check_is_EvaluationOperator(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_EvaluationOperator(v, varname)
