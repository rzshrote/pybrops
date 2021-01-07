class GenomicModelCalibrationOperator:
    """docstring for GenomicModelCalibrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(GenomicModelCalibrationOperator, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def calibrate(self, t_cur, t_max, geno, bval, **kwargs):
        """
        Calibrate genomic models using genotype and phenotype data.

        Parameters
        ----------
        geno : dict
        bval : dict

        Returns
        -------
        out : tuple
            (gmod_new, misc)
        """
        raise NotImplementedError("method is abstract")



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
