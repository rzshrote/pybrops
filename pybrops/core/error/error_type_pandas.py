from typing import Any
import pandas

################################################################################
########################## isinstance check functions ##########################
################################################################################
def check_is_pandas_df(v: Any, vname: str) -> None:
    if not isinstance(v, pandas.DataFrame):
        raise TypeError("variable '{0}' must be of type 'pandas.DataFrame'".format(vname))
