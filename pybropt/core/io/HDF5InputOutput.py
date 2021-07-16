class HDF5InputOutput:
    """docstring for HDF5InputOutput."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(HDF5InputOutput, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Matrix File I/O ####################
    def to_hdf5(self, filename, groupname):
        """
        Write object to an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name to which to write.
        groupname : str or None
            HDF5 group name under which object data is stored.
            If None, object is written to the base HDF5 group.
        """
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(cls, filename, groupname):
        """
        Read object from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which object data is stored.
            If None, object is read from base HDF5 group.

        Returns
        -------
        obj : cls
            A genotype matrix read from file.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_HDF5InputOutput(v):
    return isinstance(v, HDF5InputOutput)

def check_is_HDF5InputOutput(v, vname):
    if not isinstance(v, HDF5InputOutput):
        raise TypeError("variable '{0}' must be a HDF5InputOutput".format(vname))

def cond_check_is_HDF5InputOutput(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_HDF5InputOutput(v, vname)
