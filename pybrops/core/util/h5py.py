"""
Module containing utility functions for handling HDF5 files.
"""

import numpy
import h5py

__all__ = [
    "save_dict_to_hdf5",
]

# ruthlessly stolen/based on:
# https://codereview.stackexchange.com/questions/120802/recursively-save-python-dictionaries-to-hdf5-files-using-h5py/121308
def save_dict_to_hdf5(
        h5file: h5py.File, 
        groupname: str, 
        in_dict: dict
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
        If a field in 'in_dict' is None, skip the field; do not create a group
        for key associated with None item.
    """
    for key, item in in_dict.items():                           # for each item in dictionary
        if item is None:                                        # if item is None
            continue                                            # skip to next loop iteration
        iswritable = isinstance(                                # determine if item is writable
            item,
            (numpy.ndarray,bytes,str,int,float,numpy.floating,numpy.integer,numpy.bool_)
        )
        if iswritable:                                          # if item is writeable
            h5file.create_dataset(groupname + key, data = item) # write to dataset
        elif isinstance(item, dict):                            # else if is dictionary
            save_dict_to_hdf5(                                  # recursively call function
                h5file,                                         # input file
                groupname + key + '/',                          # append key to path
                item                                            # input dictionary
            )
        else:
            raise ValueError("Cannot save {0}: {1} type".format(key, type(item)))
