from pandas import DataFrame

################################################################################
########################## isinstance check functions ##########################
################################################################################
def check_is_pandas_DataFrame(v: object, vname: str) -> None:
    """
    Check whether a Python object is a ``pandas.DataFrame``. Raise error if not.

    Parameters
    ----------
    v : object
        Python object to check.
    vname : str
        Name associated with the Python object.
    """
    if not isinstance(v, DataFrame):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,DataFrame.__name__,type(v).__name__))
