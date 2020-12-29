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
def cond_check_is_callable(v, vname):
    generic_cond_check_hasattr(v, vname, "__call__")

def cond_check_is_iterable(v, vname):
    generic_cond_check_hasattr(v, vname, "__iter__")
