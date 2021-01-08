from . import GenomicModelCalibrationOperator

class TrueGenomicModelCalibrationOperator(GenomicModelCalibrationOperator):
    """docstring for TrueGenomicModelCalibrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(TrueGenomicModelCalibrationOperator, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def calibrate(self, t_cur, t_max, geno, bval, gmod, **kwargs):
        """
        Calibrate genomic models using genotype and phenotype data.

        Parameters
        ----------
        t_cur : int
        t_max : int
        geno : dict
        bval : dict
        gmod : dict

        Returns
        -------
        out : tuple
            (gmod_new, misc)
        """
        # shallow copy the input dictionary
        gmod_new = dict(gmod)

        # copy pointers to true genomic model to mandatory keys
        gmod_new["cand"] = gmod_new["true"]
        gmod_new["main"] = gmod_new["true"]
        gmod_new["queue"] = [gmod_new["true"] for i in range(len(gmod_new["queue"]))]

        # nothing in misc
        misc = {}

        return gmod_new, misc



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenomicModelCalibrationOperator(v):
    return isinstance(v, GenomicModelCalibrationOperator)

def check_is_GenomicModelCalibrationOperator(v, vname):
    if not isinstance(v, GenomicModelCalibrationOperator):
        raise TypeError("variable '{0}' must be a GenomicModelCalibrationOperator".format(vname))

def cond_check_is_GenomicModelCalibrationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GenomicModelCalibrationOperator(v, vname)
