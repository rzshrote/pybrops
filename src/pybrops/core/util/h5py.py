"""
Module containing utility functions for handling HDF5 files.
"""

import numpy
import h5py

from pybrops.core.error.error_type_h5py import check_is_h5py_File
from pybrops.core.error.error_value_h5py import check_h5py_File_is_writable

__all__ = [
    "save_dict_to_hdf5",
]

# writable items
writable_classes = (
    numpy.ndarray,
    bytes,
    str,
    int,
    float,
    numpy.floating,
    numpy.integer,
    numpy.bool_
)

def h5py_File_is_writable(stream: h5py.File) -> bool:
    """
    Determine if a ``h5py.File`` is writable.

    Parameters
    ----------
    stream : h5py.File
        HDF5 file stream to check.

    Returns
    -------
    out : bool
        Whether the HDF5 file is writable.
    """
    streammode = stream.file.mode
    if streammode in ("r+", "w", "w-", "x", "a"):
        return True
    return False

# ruthlessly stolen/based on:
# https://codereview.stackexchange.com/questions/120802/recursively-save-python-dictionaries-to-hdf5-files-using-h5py/121308
def save_dict_to_hdf5(
        h5file: h5py.File, 
        groupname: str, 
        in_dict: dict,
        overwrite: bool = True,
    ) -> None:
    """
    Recursively save dictionary contents to an open HDF5 file.

    Parameters
    ----------
    h5file : h5py.File
        An open HDF5 file.
    groupname : str
        String representation of group name. Must be terminated by '/'.
    in_dict : dict
        Input dictionary to save to HDF5 file.
        If a field in ``in_dict`` is ``None``, skip the field; do not create a group
        for key associated with ``None`` item.
    overwrite : bool
        Whether to overwrite fields
    """
    # type checks
    check_is_h5py_File(h5file, "h5file")
    check_h5py_File_is_writable(h5file, "h5file")

    # for each item in dictionary
    for key, item in in_dict.items():
        # if item is None, skip to the next loop iteration
        if item is None:
            continue
        
        # create field name
        fieldname = groupname + key
        
        # if item is writeable
        if isinstance(item, writable_classes):
            if (fieldname in h5file) and overwrite:
                del h5file[fieldname]
            h5file.create_dataset(fieldname, data = item)
        
        # else if is dictionary
        elif isinstance(item, dict):
            save_dict_to_hdf5(h5file, fieldname + '/', item)
        
        # else raise error
        else:
            raise ValueError("Cannot save {0}: {1} type".format(key, type(item)))
