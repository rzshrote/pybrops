"""
Error checking functions for ``h5py.File`` values.
"""

import h5py

def check_h5py_File_is_writable(stream: h5py.File, streamname: str) -> None:
    """
    Check if a ``h5py.File`` is writable, otherwise raise error.

    Parameters
    ----------
    stream : h5py.File
        HDF5 file stream to check.
    """
    streammode = stream.file.mode
    if streammode not in ("r+", "w", "w-", "x", "a"):
        raise ValueError(
            "h5py.File '{0}' must be writable but received file mode '{1}'".format(
                streamname,
                streammode
            )
        )

def check_h5py_File_has_group(stream: h5py.File, streamname: str, group: str) -> None:
    """
    Check if a ``h5py.File`` contains a required group.

    Parameters
    ----------
    stream : h5py.File
        Input ``h5py.File``.
    streamname : str
        Name of the ``h5py.File`` variable.
    group : str
        A required group that the ``h5py.File`` must contain.
    """
    if group not in stream:
        raise LookupError("h5py.File '{0}' must have group '{1}'".format(streamname,group))

def check_h5py_File_has_groups(stream: h5py.File, streamname: str, *args: tuple[str,...]) -> None:
    """
    Check if a ``h5py.File`` contains required groups.

    Parameters
    ----------
    stream : h5py.File
        Input ``h5py.File``.
    streamname : str
        Name of the ``h5py.File`` variable.
    args : tuple
        Required groups that the ``h5py.File`` must contain.
    """
    for group in args:
        if group not in stream:
            raise LookupError("h5py.File '{0}' must have group '{1}'".format(streamname,group))

