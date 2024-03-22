"""
Module containing utility functions for handling HDF5 files.
"""

import numpy
import h5py

from pybrops.core.error.error_type_h5py import check_is_h5py_File
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_value_h5py import check_h5py_File_is_readable
from pybrops.core.error.error_value_h5py import check_h5py_File_is_writable

__all__ = [
    "h5py_File_write_dict",
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

def h5py_File_is_writable(h5file: h5py.File) -> bool:
    """
    Determine if a ``h5py.File`` is writable.

    Parameters
    ----------
    h5file : h5py.File
        HDF5 file stream to check.

    Returns
    -------
    out : bool
        Whether the HDF5 file is writable.
    """
    return h5file.file.mode in ("r+", "w", "w-", "x", "a")

def h5py_File_is_readable(stream: h5py.File) -> bool:
    """
    Determine if a ``h5py.File`` is readable.

    Parameters
    ----------
    h5file : h5py.File
        HDF5 file stream to check.

    Returns
    -------
    out : bool
        Whether the HDF5 file is readable.
    """
    return stream.file.mode in ("r", "r+", "a")

def h5py_File_has_group(h5file: h5py.File, groupname: str) -> bool:
    """
    Determine if an ``h5py.File`` has a group.

    Parameters
    ----------
    h5file : h5py.File
        An HDF5 file stream to check.
    groupname : str
        Name of the group to check in ``h5file``
    
    Returns
    -------
    out : bool
        Whether ``groupname`` is in ``h5file``.
    """
    return groupname in h5file

#######################
### Write functions ###

# ruthlessly stolen/based on:
# https://codereview.stackexchange.com/questions/120802/recursively-save-python-dictionaries-to-hdf5-files-using-h5py/121308
def h5py_File_write_dict(h5file: h5py.File, groupname: str, in_dict: dict, overwrite: bool = True) -> None:
    """
    Recursively save dictionary contents to an open HDF5 file.

    Parameters
    ----------
    h5file : h5py.File
        An open, writable HDF5 file stream.
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
    check_h5py_File_is_writable(h5file)

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
            h5py_File_write_dict(h5file, fieldname + '/', item)
        
        # else raise error
        else:
            raise ValueError("Cannot save {0}: {1} type".format(key, type(item)))

######################
### Read functions ###

def h5py_File_read_int(h5file: h5py.File, fieldname: str) -> int:
    """
    Read an ``int`` from a file.

    Parameters
    ----------
    h5file : h5py.File
        An open, readable HDF5 file stream.
    fieldname : str
        Name of the field from which to read.
    
    Returns
    -------
    out : int
        A ``int`` read from the HDF5 file.
    """
    # type checks
    check_is_h5py_File(h5file, "h5file")
    check_h5py_File_is_readable(h5file)
    check_is_str(fieldname, "fieldname")

    # extract field as int
    out = int(h5file[fieldname][()])

    return out

def h5py_File_read_ndarray(h5file: h5py.File, fieldname: str) -> numpy.ndarray:
    """
    Read a ``numpy.ndarray`` from a file as is. Do not perform any type conversions.

    Bare-bones function. Does not perform any type checks on inputs.

    Parameters
    ----------
    h5file : h5py.File
        An open, readable HDF5 file stream.
    fieldname : str
        Name of the field from which to read.
    
    Returns
    -------
    out : numpy.ndarray
        An array read from the HDF5 file.
    """
    # type checks
    check_is_h5py_File(h5file, "h5file")
    check_h5py_File_is_readable(h5file)
    check_is_str(fieldname, "fieldname")

    # extract field as numpy.ndarray
    out = h5file[fieldname][()]

    return out

def h5py_File_read_ndarray_int(h5file: h5py.File, fieldname: str) -> numpy.ndarray:
    """
    Read a ``numpy.ndarray`` from a file. If the datatype is not ``int``, convert to ``int``.

    Bare-bones function. Does not perform any type checks on inputs.

    Parameters
    ----------
    h5file : h5py.File
        An open, readable HDF5 file stream.
    fieldname : str
        Name of the field from which to read.
    
    Returns
    -------
    out : numpy.ndarray
        An ``int`` array read from the HDF5 file.
    """
    # type checks
    check_is_h5py_File(h5file, "h5file")
    check_h5py_File_is_readable(h5file)
    check_is_str(fieldname, "fieldname")

    # extract field as numpy.ndarray
    out = h5file[fieldname][()]

    # if the ndarray is not integer, convert to integer
    if out.dtype != int:
        out = out.astype(int)
    
    return out

def h5py_File_read_ndarray_int8(h5file: h5py.File, fieldname: str) -> numpy.ndarray:
    """
    Read a ``numpy.ndarray`` from a file. If the datatype is not ``int8``, convert to ``int8``.

    Bare-bones function. Does not perform any type checks on inputs.

    Parameters
    ----------
    h5file : h5py.File
        An open, readable HDF5 file stream.
    fieldname : str
        Name of the field from which to read.
    
    Returns
    -------
    out : numpy.ndarray
        An ``int8`` array read from the HDF5 file.
    """
    # type checks
    check_is_h5py_File(h5file, "h5file")
    check_h5py_File_is_readable(h5file)
    check_is_str(fieldname, "fieldname")

    # extract field as numpy.ndarray
    out = h5file[fieldname][()]

    # if the ndarray is not int8, convert to int8
    if out.dtype != "int8":
        out = out.astype("int8")
    
    return out

def h5py_File_read_ndarray_utf8(h5file: h5py.File, fieldname: str) -> numpy.ndarray:
    """
    Read a ``numpy.ndarray`` from a file. If the datatype is ``bytes``, convert 
    to a utf-8 encoded ``str`` format, otherwise keep as is.

    Bare-bones function. Does not perform any type checks on inputs.

    Parameters
    ----------
    h5file : h5py.File
        An open, readable HDF5 file stream.
    fieldname : str
        Name of the field from which to read.
    
    Returns
    -------
    out : numpy.ndarray
        An ``object`` array containing strings read from the HDF5 file.
    """
    # type checks
    check_is_h5py_File(h5file, "h5file")
    check_h5py_File_is_readable(h5file)
    check_is_str(fieldname, "fieldname")

    # extract field as numpy.ndarray
    out = h5file[fieldname][()]
    
    # convert array elements to unicode / any
    out = numpy.array([s.decode("utf-8") if isinstance(s,bytes) else s for s in out], dtype = object)
    
    return out

def h5py_File_read_utf8(h5file: h5py.File, fieldname: str) -> str:
    """
    Read an ``str`` from a file.

    Parameters
    ----------
    h5file : h5py.File
        An open, readable HDF5 file stream.
    fieldname : str
        Name of the field from which to read.
    
    Returns
    -------
    out : str
        A ``str`` read from the HDF5 file.
    """
    # type checks
    check_is_h5py_File(h5file, "h5file")
    check_h5py_File_is_readable(h5file)
    check_is_str(fieldname, "fieldname")

    # extract field as bytes
    out = h5file[fieldname][()]

    # convert bytes to utf-8 string format
    out = out.decode("utf-8")

    return out

def h5py_File_read_dict(h5file: h5py.File, fieldname: str) -> dict:
    """
    Read an ``dict`` from a file.

    Parameters
    ----------
    h5file : h5py.File
        An open, readable HDF5 file stream.
    fieldname : str
        Name of the field from which to read.
    
    Returns
    -------
    out : dict
        A ``dict`` read from the HDF5 file.
    """
    # type checks
    check_is_h5py_File(h5file, "h5file")
    check_h5py_File_is_readable(h5file)
    check_is_str(fieldname, "fieldname")

    # create empty dict
    out = {}

    # get view of dataset
    view = h5file[fieldname]

    # for each field in the view, extract data
    for key in view.keys():
        out[key] = view[key][()]

    return out

