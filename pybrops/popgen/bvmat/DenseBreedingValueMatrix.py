"""
Module implementing matrix routines and associated error checking routines
for dense breeding value matrices.
"""

import copy
from numbers import Real
from typing import Optional, Union
import numpy
import h5py

from pybrops.core.error.error_type_python import check_is_array_like
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_io_h5py import check_group_in_hdf5
from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.mat.DenseTaxaTraitMatrix import DenseTaxaTraitMatrix
from pybrops.core.util.h5py import save_dict_to_hdf5
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix

class DenseBreedingValueMatrix(DenseTaxaTraitMatrix,BreedingValueMatrix):
    """
    The DenseBreedingValueMatrix class uses a dense matrix to represent a
    Multivariate Breeding Value.

    Notes
    -----
    All elements within a BreedingValueMatrix are mean-centered and scaled to
    unit variance for each trait.

    .. math::
        BV = \\frac{X - \\mu}{\\sigma}

    Where:

    - :math:`BV` is the breeding value.
    - :math:`X` is the phenotype value.
    - :math:`\\mu` is the mean (location) for :math:`X`.
    - :math:`\\sigma` is the standard deviation (scale) for :math:`X`.

    Phenotype values can be reconstituted using:

    .. math::
        X = \\sigma BV + \\mu
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            location: Union[numpy.ndarray,Real] = 0.0, 
            scale: Union[numpy.ndarray,Real] = 1.0, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        BreedingValueMatrix constructor

        Parameters
        ----------
        mat : numpy.ndarray
            An array of breeding values of shape ``(n,t)``.
            It is the responsibility of the user to ensure that the means and 
            standard deviations of this array along the ``taxa`` axis are 0 and
            1, respectively, if the breeding values are with respect to the
            individuals in the breeding value matrix.
        location : numpy.ndarray, Real
            A numpy.ndarray of shape ``(t,)`` containing breeding value locations.
            If given a Real, create a numpy.ndarray of shape ``(t,)`` filled with the provided value.
        scale : numpy.ndarray, Real
            A numpy.ndarray of shape ``(t,)`` containing breeding value scales.
            If given a Real, create a numpy.ndarray of shape ``(t,)`` filled with the provided value.
        taxa : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.
        taxa_grp : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.
        trait : numpy.ndarray, None
            A numpy.ndarray of shape ``(t,)`` containing trait names.
            If ``None``, do not store any trait name information.
        kwargs : dict
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
        out : DenseBreedingValueMatrix
            A copy of the matrix.
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
            A dictionary of objects already copied during the current copying
            pass.

        Returns
        -------
        out : DenseBreedingValueMatrix
            A deep copy of the matrix.
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
    @DenseTaxaTraitMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set raw matrix"""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim(value, "mat", 2)
        self._mat = value

    @property
    def location(self) -> numpy.ndarray:
        """Mean of the phenotype values used to calculate breeding values."""
        return self._location
    @location.setter
    def location(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set the mean of the phenotype values used to calculate breeding values"""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "location", 1)
            check_ndarray_axis_len(value, "location", 0, self.ntrait)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.ntrait)
        else:
            raise TypeError("variable 'location' must be of type 'numpy.ndarray' or 'Real'")
        self._location = value
    @location.deleter
    def location(self) -> None:
        """Delete the mean of the phenotype values used to calculate breeding values"""
        del self._location
    
    @property
    def scale(self) -> numpy.ndarray:
        """Standard deviation of the phenotype values used to calculate breeding values."""
        return self._scale
    @scale.setter
    def scale(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set the standard deviation of the phenotype values used to calculate breeding values"""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "scale", 1)
            check_ndarray_axis_len(value, "scale", 0, self.ntrait)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.ntrait)
        else:
            raise TypeError("variable 'scale' must be of type 'numpy.ndarray' or 'Real'")
        self._scale = value
    @scale.deleter
    def scale(self) -> None:
        """Delete the standard deviation of the phenotype values used to calculate breeding values"""
        del self._scale

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    # FIXME: super adjoin, delete, insert, select, ... for location, scale bug

    def select_taxa(self, indices, **kwargs: dict):
        """
        Select certain values from the Matrix along the taxa axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        kwargs : dict
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
    def targmax(self) -> numpy.ndarray:
        """
        Return indices of the maximum values for each trait column (along the taxa axis).

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing indices of maximum
            values along the taxa axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._mat.argmax(axis = self.taxa_axis)    # get argument maximum
        return out

    def targmin(self) -> numpy.ndarray:
        """
        Return indices of the minimum values for each trait column (along the taxa axis).

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing indices of minimum
            values along the taxa axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._mat.argmin(axis = self.taxa_axis)    # get argument minimum
        return out

    def tmax(self, descale: bool = False) -> numpy.ndarray:
        """
        Return the maximum for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : bool, default = False
            Whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing maximum values along the taxa
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._mat.max(axis = self.taxa_axis)   # get maximum
        if descale:
            out *= self._scale
            out += self._location
        return out

    def tmean(self, descale: bool = False) -> numpy.ndarray:
        """
        Return the mean for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : bool, default = False
            Whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing maximum values along the taxa
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._location if descale else self._mat.mean(axis = self.taxa_axis) # get mean
        return out

    def tmin(self, descale: bool = False) -> numpy.ndarray:
        """
        Return the minimum for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : bool, default = False
            Whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing minimum values along the
            taxa axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._mat.min(axis = self.taxa_axis)   # get minimum
        if descale:
            out *= self._scale
            out += self._location
        return out

    def trange(self, descale: bool = False) -> numpy.ndarray:
        """
        Return the range for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : bool, default = False
            Whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing range values along the taxa
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = numpy.ptp(self._mat, axis = self.taxa_axis)    # get range
        if descale:
            out *= self._scale
        return out

    def tstd(self, descale: bool = False) -> numpy.ndarray:
        """
        Return the standard deviation for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : bool, default = False
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing standard deviation values
            along the taxa axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._scale if descale else self._mat.std(axis = self.taxa_axis) # get standard deviation
        return out

    def tvar(self, descale: bool = False) -> numpy.ndarray:
        """
        Return the variance for each trait column (along the taxa axis).

        Parameters
        ----------
        descale : bool, default = False
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing variance values along the taxa
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._scale**2 if descale else self._mat.var(axis = self.taxa_axis) # get variance
        return out

    def descale(self) -> numpy.ndarray:
        """
        Transform values within the BreedingValueMatrix back to their de-scaled
        and de-centered values

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(n,t)`` containing de-scaled and de-centered
            values.

            Where:

            - ``n`` is the number of taxa.
            - ``t`` is the number of traits.
        """
        return (self._scale * self._mat) + self._location

    ################### Matrix File I/O ####################
    def to_hdf5(
            self, 
            filename: str, 
            groupname: Optional[str] = None
        ) -> None:
        """
        Write GenotypeMatrix to an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name to which to write.
        groupname : str, None
            HDF5 group name under which GenotypeMatrix data is stored.
            If ``None``, GenotypeMatrix is written to the base HDF5 group.
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
    def from_hdf5(
            cls, 
            filename: str, 
            groupname: Optional[str] = None
        ) -> 'DenseBreedingValueMatrix':
        """
        Read DenseBreedingValueMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which DenseBreedingValueMatrix data is stored.
            If ``None``, DenseBreedingValueMatrix is read from base HDF5 group.

        Returns
        -------
        gmat : DenseBreedingValueMatrix
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
    def from_numpy(
            cls, 
            a: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseBreedingValueMatrix':
        """
        Construct a DenseBreedingValueMatrix from a numpy.ndarray.
        Calculates mean-centering and scaling to unit variance.

        Parameters
        ----------
        a : numpy.ndarray
            A ``float64`` matrix of shape ``(n,t)``.

            Where:

            - ``n`` is the number of taxa.
            - ``t`` is the number of traits.
        taxa : numpy.ndarray
            An array of taxa names.
        taxa_grp : numpy.ndarray
            An array of taxa groups.
        trait : numpy.ndarray
            An array of trait names.

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
def is_DenseBreedingValueMatrix(v: object) -> bool:
    """
    Determine whether an object is a DenseBreedingValueMatrix.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseBreedingValueMatrix object instance.
    """
    return isinstance(v, DenseBreedingValueMatrix)

def check_is_DenseBreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseBreedingValueMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseBreedingValueMatrix):
        raise TypeError("variable '{0}' must be a DenseBreedingValueMatrix".format(vname))
