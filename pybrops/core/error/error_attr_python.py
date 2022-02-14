from . import generic_check_hasattr

from . import generic_default_cond
from . import generic_cond_check_hasattr

### read/write ###
def error_readonly(vname):
    raise AttributeError("variable '{0}' is read-only".format(vname))

################################################################################
########################### attribute check functions ##########################
################################################################################
# do not use callable() due to removal in Python 3.0
def check_is_callable(v, vname):
    generic_check_hasattr(v, vname, "__call__")

def check_is_iterable(v, vname):
    generic_check_hasattr(v, vname, "__iter__")

################################################################################
#################### conditional attribute check functions #####################
################################################################################
def cond_check_is_callable(v, vname, cond = generic_default_cond):
    generic_cond_check_hasattr(v, vname, "__call__", cond)

def cond_check_is_iterable(v, vname, cond = generic_default_cond):
    generic_cond_check_hasattr(v, vname, "__iter__", cond)
