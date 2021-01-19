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
    def evaluate(self, t_cur, t_max, pgvmat, gmod_true, **kwargs):
        """
        Parameters
        ----------
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        pgvmat : PhasedGenotypeVariantMatrix
            Genotypes to evaluate.
        gmod_true : GenomicModel
            True genomic model.

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
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_EvaluationOperator(v):
    return isinstance(v, EvaluationOperator)

def check_is_EvaluationOperator(v, vname):
    if not isinstance(v, EvaluationOperator):
        raise TypeError("variable '{0}' must be a EvaluationOperator".format(vname))

def cond_check_is_EvaluationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_EvaluationOperator(v, vname)
