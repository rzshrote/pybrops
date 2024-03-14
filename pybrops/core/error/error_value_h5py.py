"""
Error checking functions for ``h5py.File`` values.
"""

import h5py

def check_h5py_File_is_writable(h5file: h5py.File) -> None:
    """
    Check if a ``h5py.File`` is writable, otherwise raise error.

    Parameters
    ----------
    h5file : h5py.File
        HDF5 file stream to check.
    """
    if h5file.file.mode not in ("r+", "w", "w-", "x", "a"):
        raise ValueError(
            "``h5py.File`` (file: {0}) must be writable but received file mode ``{1}``".format(
                h5file.filename,
                h5file.file.mode,
            )
        )

def check_h5py_File_is_readable(h5file: h5py.File) -> None:
    """
    Check if a ``h5py.File`` is readable, otherwise raise error.

    Parameters
    ----------
    h5file : h5py.File
        HDF5 file stream to check.
    streamname : str
        Name of the ``h5py.File`` variable.
    """
    if h5file.file.mode not in ("r", "r+", "a"):
        raise ValueError(
            "``h5py.File`` (file: {0}) must be readable but received file mode ``{1}``".format(
                h5file.filename,
                h5file.file.mode,
            )
        )

def check_h5py_File_has_group(h5file: h5py.File, group: str) -> None:
    """
    Check if a ``h5py.File`` contains a required group.

    Parameters
    ----------
    h5file : h5py.File
        Input ``h5py.File``.
    streamname : str
        Name of the ``h5py.File`` variable.
    group : str
        A required group that the ``h5py.File`` must contain.
    """
    if group not in h5file:
        raise LookupError(
            "``h5py.File`` (file: {0}) must have group ``{1}``".format(
                h5file.filename,
                group
            )
        )

def check_h5py_File_has_groups(h5file: h5py.File, *groups: tuple[str,...]) -> None:
    """
    Check if a ``h5py.File`` contains required groups.

    Parameters
    ----------
    h5file : h5py.File
        Input ``h5py.File``.
    streamname : str
        Name of the ``h5py.File`` variable.
    groups : tuple
        Required groups that the ``h5py.File`` must contain.
    """
    for group in groups:
        if group not in h5file:
            raise LookupError(
                "``h5py.File`` (file: {0}) must have group ``{1}``".format(
                    h5file.filename,
                    group
                )
            )

