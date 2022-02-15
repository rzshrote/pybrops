from pybrops.core.mat import TaxaMatrix

class BreedingValueMatrix(TaxaMatrix):
    """docstring for BreedingValueMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        BreedingValueMatrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(BreedingValueMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Breeding Value Data ##################
    def trait():
        doc = "The trait property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    trait = property(**trait())

    def ntrait():
        doc = "The ntrait property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    ntrait = property(**ntrait())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Matrix File I/O ####################
    @staticmethod
    def from_hdf5(filename, groupname):
        """
        Read GenotypeMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is read from base HDF5 group.

        Returns
        -------
        gmat : GenotypeMatrix
            A genotype matrix read from file.
        """
        raise NotImplementedError("method is abstract")

    def to_hdf5(self, filename, groupname):
        """
        Write GenotypeMatrix to an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name to which to write.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is written to the base HDF5 group.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingValueMatrix(v):
    return isinstance(v, BreedingValueMatrix)

def check_is_BreedingValueMatrix(v, vname):
    if not isinstance(v, BreedingValueMatrix):
        raise TypeError("variable '{0}' must be a BreedingValueMatrix".format(vname))

def cond_check_is_BreedingValueMatrix(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_BreedingValueMatrix(v, vname)
