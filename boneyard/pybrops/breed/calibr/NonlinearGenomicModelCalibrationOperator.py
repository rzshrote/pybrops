from . import GenomicModelCalibrationOperator

class NonlinearGenomicModelCalibrationOperator(GenomicModelCalibrationOperator):
    """docstring for NonlinearGenomicModelCalibrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict):
        super(NonlinearGenomicModelCalibrationOperator, self).__init__(
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_NonlinearGenomicModelCalibrationOperator(v):
    return isinstance(v, NonlinearGenomicModelCalibrationOperator)

def check_is_NonlinearGenomicModelCalibrationOperator(v, varname):
    if not isinstance(v, NonlinearGenomicModelCalibrationOperator):
        raise TypeError("'%s' must be a NonlinearGenomicModelCalibrationOperator." % varname)

def cond_check_is_NonlinearGenomicModelCalibrationOperator(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_NonlinearGenomicModelCalibrationOperator(v, varname)
