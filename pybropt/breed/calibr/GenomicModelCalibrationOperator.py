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
    def calibrate(self, **kwargs):
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenomicModelCalibrationOperator(v):
    return isinstance(v, GenomicModelCalibrationOperator)

def check_is_GenomicModelCalibrationOperator(v, varname):
    if not isinstance(v, GenomicModelCalibrationOperator):
        raise TypeError("'%s' must be a GenomicModelCalibrationOperator." % varname)

def cond_check_is_GenomicModelCalibrationOperator(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GenomicModelCalibrationOperator(v, varname)
