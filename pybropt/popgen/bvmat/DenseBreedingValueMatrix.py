import numpy

from pybropt.core.mat import DenseTaxaTraitMatrix
from . import BreedingValueMatrix
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_ndim
from pybropt.core.util import save_dict_to_hdf5

class DenseBreedingValueMatrix(DenseTaxaTraitMatrix,BreedingValueMatrix):
    """Dense breeding value matrix implementation."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, trait = None, **kwargs):
        """
        BreedingValueMatrix constructor

        Parameters
        ----------
        mat : numpy.ndarray
            A float64 matrix of breeding values of shape (n, t).
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.

        """
        super(DenseBreedingValueMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Breeding Value Data ##################
    def mat():
        doc = "Raw matrix property."
        def fget(self):
            """Get raw matrix"""
            return self._mat
        def fset(self, value):
            """Set raw matrix"""
            check_is_ndarray(value, "mat")
            check_ndarray_ndim(value, "mat", 2)
            self._mat = value
        def fdel(self):
            """Delete raw matrix"""
            del self._mat
        return locals()
    mat = property(**mat())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############## Matrix summary statistics ###############
    def targmax(self):
        """
        Return indices of the maximum values along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape (t,) containing indices of maximum values
            along the trait axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.argmax(axis = 0)    # get argument maximum
        return out

    def targmin(self):
        """
        Return indices of the minimum values along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape (t,) containing indices of minimum values
            along the trait axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.argmin(axis = 0)    # get argument minimum
        return out

    def tmax(self):
        """
        Return the maximum along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing maximum values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.max(axis = 0)   # get maximum
        return out

    def tmean(self):
        """
        Return the mean along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing maximum values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.mean(axis = 0)  # get mean
        return out

    def tmin(self):
        """
        Return the minimum along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape (t,) containing minimum values along the
            trait axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.min(axis = 0)   # get minimum
        return out

    def trange(self):
        """
        Return the range along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing variance values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        out = numpy.ptp(self._mat, axis = 0)    # get range
        return out

    def tstd(self):
        """
        Return the standard deviation along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing standard deviation values along
            the trait axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.std(axis = 0)   # get standard deviation
        return out

    def tvar(self):
        """
        Return the variance along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing variance values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.var(axis = 0)   # get variance
        return out

    ################### Matrix File I/O ####################
    def to_hdf5(self, filename, groupname = None):
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
        h5file = h5py.File(filename, "a")                       # open HDF5 in write mode
        ######################################################### process groupname argument
        if isinstance(groupname, str):                          # if we have a string
            if groupname[-1] != '/':                            # if last character in string is not '/'
                groupname += '/'                                # add '/' to end of string
        elif groupname is None:                                 # else if groupname is None
            groupname = ""                                      # empty string
        else:                                                   # else raise error
            raise TypeError("'groupname' must be of type str or None")
        ######################################################### populate HDF5 file
        data_dict = {                                           # data dictionary
            "mat": self.mat,
            "taxa": self.taxa,
            "taxa_grp": self.taxa_grp,
            "trait": self.trait
        }
        save_dict_to_hdf5(h5file, groupname, data_dict)         # save data
        ######################################################### write conclusion
        h5file.close()                                          # close the file

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(cls, filename, groupname = None):
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
        check_file_exists(filename)                             # check file exists
        h5file = h5py.File(filename, "r")                       # open HDF5 in read only
        ######################################################### process groupname argument
        if isinstance(groupname, str):                          # if we have a string
            check_group_in_hdf5(groupname, h5file, filename)    # check that group exists
            if groupname[-1] != '/':                            # if last character in string is not '/'
                groupname += '/'                                # add '/' to end of string
        elif groupname is None:                                 # else if groupname is None
            groupname = ""                                      # empty string
        else:                                                   # else raise error
            raise TypeError("'groupname' must be of type str or None")
        ######################################################### check that we have all required fields
        required_fields = ["mat"]                               # all required arguments
        for field in required_fields:                           # for each required field
            fieldname = groupname + field                       # concatenate base groupname and field
            check_group_in_hdf5(fieldname, h5file, filename)    # check that group exists
        ######################################################### read data
        data_dict = {                                           # output dictionary
            "mat": None,
            "taxa": None,
            "taxa_grp": None,
            "trait": None
        }
        for field in data_dict.keys():                          # for each field
            fieldname = groupname + field                       # concatenate base groupname and field
            if fieldname in h5file:                             # if the field exists in the HDF5 file
                data_dict[field] = h5file[fieldname][:]         # read array
        ######################################################### read conclusion
        h5file.close()                                          # close file
        data_dict["taxa"] = numpy.object_(                      # convert taxa strings from byte to utf-8
            [s.decode("utf-8") for s in data_dict["taxa"]]
        )
        data_dict["trait"] = numpy.object_(                     # convert trait string from byte to utf-8
            [s.decode("utf-8") for s in data_dict["trait"]]
        )
        ######################################################### create object
        gmat = cls(**data_dict)                                 # create object from read data
        return gmat



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseBreedingValueMatrix(v):
    return isinstance(v, DenseBreedingValueMatrix)

def check_is_DenseBreedingValueMatrix(v, varname):
    if not isinstance(v, DenseBreedingValueMatrix):
        raise TypeError("'%s' must be a DenseBreedingValueMatrix." % varname)

def cond_check_is_DenseBreedingValueMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DenseBreedingValueMatrix(v, varname)
