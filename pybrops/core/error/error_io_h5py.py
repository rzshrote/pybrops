def check_group_in_hdf5(groupname, h5file, h5filename):
    if groupname not in h5file:
        raise LookupError("{0} group does not exist in {1}".format(groupname, h5filename))
