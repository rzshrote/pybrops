import copy
import numpy
import h5py

from pybropt.core.mat import DenseTaxaTraitMatrix
from . import BreedingValueMatrix
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_ndim
from pybropt.core.util import save_dict_to_hdf5
from pybropt.core.error import check_ndarray_std_is_approx
from pybropt.core.error import check_ndarray_mean_is_approx
from pybropt.core.error import check_ndarray_axis_len
from pybropt.core.error import check_is_array_like

class DenseBreedingValueMatrix(DenseTaxaTraitMatrix,BreedingValueMatrix):
    """Dense breeding value matrix implementation."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, location, scale, taxa = None, taxa_grp = None, trait = None, **kwargs):
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
        # set location and scale parameters
        self.location = location
        self.scale = scale

    #################### Matrix copying ####################
    def __copy__(self):
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : Matrix
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            location = copy.copy(self.location),
            scale = copy.copy(self.scale),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
            trait = copy.copy(self.trait)
        )
        # copy taxa metadata
        out.taxa_grp_name = copy.copy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.copy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.copy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.copy(self.taxa_grp_len)

        return out

    def __deepcopy__(self, memo):
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict

        Returns
        -------
        out : Matrix
        """
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            location = copy.deepcopy(self.location),
            scale = copy.deepcopy(self.scale),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo),
            trait = copy.deepcopy(self.trait, memo)
        )

        # copy taxa metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name, memo)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix, memo)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix, memo)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len, memo)

        return out

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
            check_ndarray_mean_is_approx(value, "mat", 0.0, self.taxa_axis)
            check_ndarray_std_is_approx(value, "mat", 1.0, self.taxa_axis)
            self._mat = value
        def fdel(self):
            """Delete raw matrix"""
            del self._mat
        return locals()
    mat = property(**mat())

    def location():
        doc = "Mean of the phenotype values used to calculate breeding values"
        def fget(self):
            """Get the mean of the phenotype values used to calculate breeding values"""
            return self._location
        def fset(self, value):
            """Set the mean of the phenotype values used to calculate breeding values"""
            check_is_ndarray(value, "location")
            check_ndarray_ndim(value, "location", 1)
            check_ndarray_axis_len(value, "location", 0, self.ntrait)
            self._location = value
        def fdel(self):
            """Delete the mean of the phenotype values used to calculate breeding values"""
            del self._location
        return locals()
    location = property(**location())

    def scale():
        doc = "Standard deviation of the phenotype values used to calculate breeding values"
        def fget(self):
            """Get the standard deviation of the phenotype values used to calculate breeding values"""
            return self._scale
        def fset(self, value):
            """Set the standard deviation of the phenotype values used to calculate breeding values"""
            check_is_ndarray(value, "scale")
            check_ndarray_ndim(value, "scale", 1)
            check_ndarray_axis_len(value, "scale", 0, self.ntrait)
            self._scale = value
        def fdel(self):
            """Delete the standard deviation of the phenotype values used to calculate breeding values"""
            del self._scale
        return locals()
    scale = property(**scale())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    # FIXME: super adjoin, delete, insert, select, ... for location, scale bug

    def select_taxa(self, indices, **kwargs):
        """
        Select certain values from the Matrix along the taxa axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            The output Matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        # check for array_like
        check_is_array_like(indices, "indices")

        # get values
        mat = self.descale()        # get descaled values
        taxa = self._taxa
        taxa_grp = self._taxa_grp

        # select values
        mat = numpy.take(mat, indices, axis = self.taxa_axis)
        if taxa is not None:
            taxa = numpy.take(taxa, indices, axis = 0)
        if taxa_grp is not None:
            taxa_grp = numpy.take(taxa_grp, indices, axis = 0)

        # re-calculate breeding values
        location = mat.mean(0)          # recalculate location
        scale = mat.std(0)              # recalculate scale
        mat = (mat - location) / scale  # mean center and scale values

        # construct output
        out = self.__class__(
            mat = mat,
            location = location,
            scale = scale,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

        return out

    ############## Matrix summary statistics ###############
    def targmax(self):
        """
        Return indices of the maximum values for each trait column (along the taxa axis).

        Returns
        -------
        out : numpy.ndarray
            An index array of shape (t,) containing indices of maximum values
            along the taxa axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.argmax(axis = self.taxa_axis)    # get argument maximum
        return out

    def targmin(self):
        """
        Return indices of the minimum values for each trait column (along the taxa axis).

        Returns
        -------
        out : numpy.ndarray
            An index array of shape (t,) containing indices of minimum values
            along the taxa axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.argmin(axis = self.taxa_axis)    # get argument minimum
        return out

    def tmax(self, descale = False):
        """
        Return the maximum for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : boolean, default = False
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing maximum values along the taxa
            axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.max(axis = self.taxa_axis)   # get maximum
        if descale:
            out *= self._scale
            out += self._location
        return out

    def tmean(self, descale = False):
        """
        Return the mean for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : boolean, default = False
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing maximum values along the taxa
            axis.
            Where:
                't' is the number of traits.
        """
        out = self._location if descale else self._mat.mean(axis = self.taxa_axis) # get mean
        return out

    def tmin(self, descale = False):
        """
        Return the minimum for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : boolean, default = False
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape (t,) containing minimum values along the
            taxa axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.min(axis = self.taxa_axis)   # get minimum
        if descale:
            out *= self._scale
            out += self._location
        return out

    def trange(self, descale = False):
        """
        Return the range for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : boolean, default = False
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing range values along the taxa
            axis.
            Where:
                't' is the number of traits.
        """
        out = numpy.ptp(self._mat, axis = self.taxa_axis)    # get range
        if descale:
            out *= self._scale
            out += self._location
        return out

    def tstd(self, descale = False):
        """
        Return the standard deviation for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : boolean, default = False
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing standard deviation values along
            the taxa axis.
            Where:
                't' is the number of traits.
        """
        out = self._scale if descale else self._mat.std(axis = self.taxa_axis) # get standard deviation
        return out

    def tvar(self, descale = False):
        """
        Return the variance for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : boolean, default = False
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing variance values along the taxa
            axis.
            Where:
                't' is the number of traits.
        """
        out = self._scale**2 if descale else self._mat.var(axis = self.taxa_axis) # get variance
        return out

    def descale(self):
        """
        Transform values within the BreedingValueMatrix back to their de-scaled
        and de-centered values

        Returns
        -------
        out : numpy.ndarray
            An array of shape (n,t) containing de-scaled and de-centered values.
            Where:
                'n' is the number of taxa.
                't' is the number of traits.
        """
        return (self._scale * self._mat) + self._location

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
            "location": self.location,
            "scale": self.scale,
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
        required_fields = ["mat", "location", "scale"]          # all required arguments
        for field in required_fields:                           # for each required field
            fieldname = groupname + field                       # concatenate base groupname and field
            check_group_in_hdf5(fieldname, h5file, filename)    # check that group exists
        ######################################################### read data
        data_dict = {                                           # output dictionary
            "mat": None,
            "location": None,
            "scale": None,
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

    @classmethod
    def from_numpy(cls, a, taxa = None, taxa_grp = None, trait = None, **kwargs):
        """
        Construct a DenseBreedingValueMatrix from a numpy.ndarray.
        Calculates mean-centering and scaling to unit variance.

        Parameters
        ----------
        a : numpy.ndarray
            A float64 matrix of shape (n,t).
            Where:
                'n' is the number of taxa.
                't' is the number of traits.
        taxa : numpy.ndarray
        taxa_grp : numpy.ndarray
        trait : numpy.ndarray

        Returns
        -------
        out : DenseBreedingValueMatrix
            Output breeding value matrix.
        """
        # check inputs
        check_ndarray_ndim(a, "a", 2)

        # calculate location and scale parameters
        location = a.mean(0)
        scale = a.std(0)

        # mean center and scale values
        mat = (a - location) / scale

        # construct output
        out = cls(
            mat = mat,
            location = location,
            scale = scale,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

        return out



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseBreedingValueMatrix(v):
    """Return whether an object is a DenseBreedingValueMatrix or not"""
    return isinstance(v, DenseBreedingValueMatrix)

def check_is_DenseBreedingValueMatrix(v, varname):
    """Raise TypeError if object is not a DenseBreedingValueMatrix"""
    if not isinstance(v, DenseBreedingValueMatrix):
        raise TypeError("'%s' must be a DenseBreedingValueMatrix." % varname)

def cond_check_is_DenseBreedingValueMatrix(v, varname, cond=(lambda s: s is not None)):
    """If object is not None, raise TypeError if object is not a DenseBreedingValueMatrix"""
    if cond(v):
        check_is_DenseBreedingValueMatrix(v, varname)
