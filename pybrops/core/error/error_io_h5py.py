"""
Module containing error subroutines related to HDF5 files.
"""

__all__ = [
    "check_group_in_hdf5"
]

def check_group_in_hdf5(groupname: str, h5file, h5filename: str):
    """
    Subroutine to check whether a given group is in an HDF5 file.
    If the group is not in the HDF5 file, raise an LookupError with a custom
    error message.

    Parameters
    ----------
    groupname : str
        Name of the group to search for in the HDF5 file.
    h5file : File
        File pointer to an HDF5 file.
    h5filename : str
        String representing a name given to ``h5file``. This does not
        necessarily need to be a valid file path, it could just be a descriptor
        for the HDF5 file (e.g. "genotype HDF5 file")
    """
    if groupname not in h5file:
        raise LookupError(
            "{0} group does not exist in {1}".format(groupname, h5filename)
        )
