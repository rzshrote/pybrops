import pandas

from . import generic_check_isinstance
from . import generic_cond_check_isinstance
from . import generic_default_cond

################################################################################
########################## isinstance check functions ##########################
################################################################################
def check_is_pandas_df(v, vname):
    generic_check_isinstance(v, vname, pandas.DataFrame)

################################################################################
#################### conditional isinstance check functions ####################
################################################################################
def cond_check_is_pandas_df(v, vname, cond = generic_default_cond):
    generic_cond_check_isinstance(v, vname, pandas.DataFrame, cond)
