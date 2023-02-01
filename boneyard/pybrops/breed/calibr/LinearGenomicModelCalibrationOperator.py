from . import GenomicModelCalibrationOperator

class LinearGenomicModelCalibrationOperator(GenomicModelCalibrationOperator):
    """docstring for LinearGenomicModelCalibrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict):
        super(LinearGenomicModelCalibrationOperator, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_LinearGenomicModelCalibrationOperator(v):
    return isinstance(v, LinearGenomicModelCalibrationOperator)

def check_is_LinearGenomicModelCalibrationOperator(v, varname):
    if not isinstance(v, LinearGenomicModelCalibrationOperator):
        raise TypeError("'%s' must be a LinearGenomicModelCalibrationOperator." % varname)

def cond_check_is_LinearGenomicModelCalibrationOperator(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_LinearGenomicModelCalibrationOperator(v, varname)
