"""
Error checking routines for h5py.
"""

import h5py

def check_is_h5py_File(v: object, vname: str) -> None:
    """
    Check whether a Python object is a ``h5py.File``.

    Parameters
    ----------
    v : object
        Python object to test.
    vname : str
        Name assigned to the Python object for an error message.
    """
    if not isinstance(v, h5py.File):
        raise TypeError(
            "variable '{0}' must be of type '{1}' but received type '{2}'".format(
                vname,
                "h5py.File",
                type(v).__name__
            )
        )
