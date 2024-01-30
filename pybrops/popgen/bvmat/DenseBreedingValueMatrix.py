"""
Module implementing matrix routines and associated error checking routines
for dense breeding value matrices.
"""

__all__ = [
    "DenseBreedingValueMatrix",
    "check_is_DenseBreedingValueMatrix",
]

import copy
from numbers import Integral, Real
from typing import Optional, Sequence, Union
import numpy
import h5py
import pandas
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_type_python import check_is_array_like
from pybrops.core.error.error_type_python import check_is_bool
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_type_python import check_is_str_or_Integral
from pybrops.core.error.error_type_python import check_is_str_or_Sequence
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group
from pybrops.core.error.error_value_numpy import check_ndarray_all_gteq
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_value_python import check_is_gteq
from pybrops.core.error.error_value_python import check_len
from pybrops.core.error.error_value_python import check_str_value
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

    ########################## Special Object Methods ##########################
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
            A ``numpy.ndarray`` of shape ``(t,)`` containing breeding value 
            locations. If given a ``Real``, create a ``numpy.ndarray`` of shape 
            ``(t,)`` filled with the provided value.
        
        scale : numpy.ndarray, Real
            A ``numpy.ndarray`` of shape ``(t,)`` containing breeding value 
            scales. If given a ``Real``, create a ``numpy.ndarray`` of shape 
            ``(t,)`` filled with the provided value.
        
        taxa : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.
        
        taxa_grp : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.
        
        trait : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(t,)`` containing trait names.
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
    def __copy__(
            self
        ) -> 'DenseBreedingValueMatrix':
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

    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'DenseBreedingValueMatrix':
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

    ############################ Object Properties #############################

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
            check_ndarray_all_gteq(value, "scale", 0)
        elif isinstance(value, Real):
            check_is_gteq(value, "scale", 0)
            value = numpy.repeat(value, self.ntrait)
        else:
            raise TypeError("variable 'scale' must be of type 'numpy.ndarray' or 'Real'")
        self._scale = value

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseBreedingValueMatrix':
        """
        Make a shallow copy of the Matrix.

        Returns
        -------
        out : DenseMatrix
            A shallow copy of the original DenseMatrix.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseBreedingValueMatrix':
        """
        Make a deep copy of the Matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseMatrix
            A deep copy of the original DenseMatrix.
        """
        return copy.deepcopy(self, memo)

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
        mat = self.unscale()        # get unscaled values
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

    def tmax(self, unscale: bool = False) -> numpy.ndarray:
        """
        Return the maximum for each trait column (along the taxa axis).

        Parameters
        ----------
        unscale : bool, default = False
            Whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing maximum values along the taxa
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._mat.max(axis = self.taxa_axis)   # get maximum
        if unscale:
            out *= self._scale
            out += self._location
        return out

    def tmean(self, unscale: bool = False) -> numpy.ndarray:
        """
        Return the mean for each trait column (along the taxa axis).

        Parameters
        ----------
        unscale : bool, default = False
            Whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing maximum values along the taxa
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._location if unscale else self._mat.mean(axis = self.taxa_axis) # get mean
        return out

    def tmin(self, unscale: bool = False) -> numpy.ndarray:
        """
        Return the minimum for each trait column (along the taxa axis).

        Parameters
        ----------
        unscale : bool, default = False
            Whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing minimum values along the
            taxa axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._mat.min(axis = self.taxa_axis)   # get minimum
        if unscale:
            out *= self._scale
            out += self._location
        return out

    def trange(self, unscale: bool = False) -> numpy.ndarray:
        """
        Return the range for each trait column (along the taxa axis).

        Parameters
        ----------
        unscale : bool, default = False
            Whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing range values along the taxa
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = numpy.ptp(self._mat, axis = self.taxa_axis)    # get range
        if unscale:
            out *= self._scale
        return out

    def tstd(self, unscale: bool = False) -> numpy.ndarray:
        """
        Return the standard deviation for each trait column (along the taxa axis).

        Parameters
        ----------
        unscale : bool, default = False
            whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing standard deviation values
            along the taxa axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._scale if unscale else self._mat.std(axis = self.taxa_axis) # get standard deviation
        return out

    def tvar(self, unscale: bool = False) -> numpy.ndarray:
        """
        Return the variance for each trait column (along the taxa axis).

        Parameters
        ----------
        unscale : bool, default = False
            whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing variance values along the taxa
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        out = self._scale**2 if unscale else self._mat.var(axis = self.taxa_axis) # get variance
        return out

    def unscale(self) -> numpy.ndarray:
        """
        Transform values within the BreedingValueMatrix back to their unscaled
        and de-centered values

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(n,t)`` containing unscaled and de-centered
            values.

            Where:

            - ``n`` is the number of taxa.
            - ``t`` is the number of traits.
        """
        return (self._scale * self._mat) + self._location

    ################### Matrix File I/O ####################
    def to_pandas(
            self, 
            taxa_col: Optional[str] = "taxa",
            taxa_grp_col: Optional[str] = "taxa_grp",
            trait_cols: Optional[Union[str,Sequence]] = "all",
            unscale: bool = False,
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Export a DenseBreedingValueMatrix to a pandas.DataFrame.

        Parameters
        ----------
        taxa_col : str, None, default = "taxa"
            Name of the column to which to write taxa names.
            If ``str``, the column is given the name in ``taxa_col``.
            If ``None``, the column is not exported.
        
        taxa_grp_col : str, None, default = "taxa_grp"
            Name of the column to which to write taxa group names.
            If ``str``, the column is given the name in ``taxa_grp_col``.
            If ``None``, the column is not exported.

        trait_cols : Sequence, str, None, default = "trait"
            Names of the trait columns to which to write breeding values.
            If ``Sequence``, column names are given by the strings in the 
            ``trait_cols`` Sequence.
            If ``str``, must be equal to ``"all"``. Use all trait names given 
            in the ``trait`` property.
            If ``None``, use numeric trait column names.
        
        unscale : bool, default = False
            whether to transform breeding values to their unscaled values.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            pandas.DataFrame.
        
        Returns
        -------
        out : pandas.DataFrame
            An output dataframe.
        """
        # type checks
        if taxa_col is not None:
            check_is_str(taxa_col, "taxa_col")
        if taxa_grp_col is not None:
            check_is_str(taxa_grp_col, "taxa_grp_col")
        if trait_cols is not None:
            if isinstance(trait_cols, str):
                check_str_value(trait_cols, "trait_cols", "all")
            elif isinstance(trait_cols, Sequence):
                check_len(trait_cols, "trait_cols", self.ntrait)
            else:
                check_is_str_or_Sequence(trait_cols, "trait_cols")
        check_is_bool(unscale, "unscale")

        # construct dictionary for labels and data
        data_dict = {}

        # process taxa_col
        if taxa_col is not None:
            data_dict[taxa_col] = self.taxa
        
        # process taxa_grp_col
        if taxa_grp_col is not None:
            data_dict[taxa_grp_col] = self.taxa_grp
        
        # process trait_cols
        if trait_cols is None:
            trait_cols = numpy.arange(self.ntrait)
        elif isinstance(trait_cols, str):
            trait_cols = numpy.arange(self.ntrait) if self.trait is None else self.trait
        
        # extract breeding values
        bv = self.unscale() if unscale else self.mat
        for i,trait in zip(range(self.ntrait),trait_cols):
            data_dict[trait] = bv[:,i]
        
        # create dataframe
        out = pandas.DataFrame(data_dict)

        return out

    def to_csv(
            self,
            filename: str,
            taxa_col: Optional[str] = "taxa",
            taxa_grp_col: Optional[str] = "taxa_grp",
            trait_cols: Optional[Union[str,Sequence]] = "trait",
            unscale: bool = False,
            sep: str = ',', 
            header: bool = True, 
            index: bool = False, 
            **kwargs: dict
        ) -> None:
        """
        Write a DenseBreedingValueMatrix to a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name to which to write.
        
        taxa_col : str, None, default = "taxa"
            Name of the column to which to write taxa names.
            If ``str``, the column is given the name in ``taxa_col``.
            If ``None``, the column is not exported.
        
        taxa_grp_col : str, None, default = "taxa_grp"
            Name of the column to which to write taxa group names.
            If ``str``, the column is given the name in ``taxa_grp_col``.
            If ``None``, the column is not exported.

        trait_cols : Sequence, str, None, default = "trait"
            Names of the trait columns to which to write breeding values.
            If ``Sequence``, column names are given by the strings in the 
            ``trait_cols`` Sequence.
            If ``str``, must be equal to ``"trait"``. Use trait names given in 
            the ``trait`` property.
            If ``None``, use numeric trait column names.
        
        unscale : bool, default = False
            whether to transform breeding values to their unscaled values.
        
        sep : str, default = ","
            Separator to use in the exported CSV file.
        
        header : bool, default = True
            Whether to save header names.
        
        index : bool, default = False
            Whether to save a row index in the exported CSV file.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a CSV.
        """
        # convert DenseBreedingValueMatrix to pandas.DataFrame
        df = self.to_pandas(
            taxa_col = taxa_col,
            taxa_grp_col = taxa_grp_col,
            trait_cols = trait_cols,
            unscale = unscale,
        )

        # export using pandas
        df.to_csv(
            path_or_buf = filename,
            sep = sep,
            header = header,
            index = index,
            **kwargs
        )

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

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    def from_pandas(
            cls, 
            df: pandas.DataFrame,
            location: Union[numpy.ndarray,Real] = 0.0, 
            scale: Union[numpy.ndarray,Real] = 1.0, 
            taxa_col: Optional[Union[str,Integral]] = "taxa",
            taxa_grp_col: Optional[Union[str,Integral]] = "taxa_grp",
            trait_cols: Optional[Union[str,Sequence]] = "infer",
            **kwargs: dict
        ) -> 'DenseBreedingValueMatrix':
        """
        Read a DenseBreedingValueMatrix from a pandas.DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.

        location : numpy.ndarray, Real, default = 0.0
            A ``numpy.ndarray`` of shape ``(t,)`` containing breeding value 
            locations. If given a ``Real``, create a ``numpy.ndarray`` of shape 
            ``(t,)`` filled with the provided value.

        scale : numpy.ndarray, Real, default = 1.0
            A ``numpy.ndarray`` of shape ``(t,)`` containing breeding value 
            scales. If given a ``Real``, create a ``numpy.ndarray`` of shape 
            ``(t,)`` filled with the provided value.

        taxa_col : str, Integral, None, default = "taxa"
            Name of the column from which to read taxa names.
            If of type ``str``, taxa names are read from the column named 
            defined by ``taxa_col``.
            If of type ``Integral``, taxa names are read from the column 
            number defined by ``taxa_col``.
            If ``None``, taxa names are not imported.
        
        taxa_grp_col : str, None, default = "taxa_grp"
            Name of the column to which to read taxa group names.
            If of type ``str``, taxa group names are read from the column named 
            defined by ``taxa_col``.
            If of type ``Integral``, taxa group names are read from the column 
            number defined by ``taxa_col``.
            If ``None``, taxa group names are not imported.

        trait_cols : Sequence, str, None, default = "trait"
            Names of the trait columns to which to read breeding values.
            If ``Sequence``, column names are given by the strings or integers 
            in the ``trait_cols`` Sequence.
            If ``str``, must be equal to ``"infer"``. Use remaining columns in 
            the input dataframe to load trait breeding values.
            If ``None``, do not load any trait breeding values.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            pandas.DataFrame.

        Returns
        -------
        out : DenseBreedingValueMatrix
            A DenseBreedingValueMatrix read from a pandas.DataFrame.
        """
        # type checks
        check_is_pandas_DataFrame(df, "df")
        if taxa_col is not None:
            check_is_str_or_Integral(taxa_col, "taxa_col")
        if taxa_grp_col is not None:
            check_is_str_or_Integral(taxa_grp_col, "taxa_grp_col")
        if trait_cols is not None:
            if isinstance(trait_cols, str):
                check_str_value(trait_cols, "trait_cols", "infer")
            elif isinstance(trait_cols, Sequence):
                pass
            else:
                check_is_str_or_Sequence(trait_cols, "trait_cols")

        ### extract data from dataframe
        colmask = numpy.full(len(df.columns), True, dtype = bool)

        # extract taxa data
        taxa = None
        if taxa_col is not None:
            taxaix = df.columns.get_loc(taxa_col) if isinstance(taxa_col, str) else taxa_col
            taxa = df.iloc[:,taxaix].to_numpy(dtype = object)
            colmask[taxaix] = False
        
        # extract taxa group data
        taxa_grp = None
        if taxa_grp_col is not None:
            taxagrpix = df.columns.get_loc(taxa_grp_col) if isinstance(taxa_grp_col, str) else taxa_grp_col
            taxa_grp = df.iloc[:,taxagrpix].to_numpy(dtype = int)
            colmask[taxagrpix] = False
        
        # for non-string Sequence, re-construct column mask
        if isinstance(trait_cols, Sequence) and not isinstance(trait_cols, str):
            colmask[:] = False
            for trait_col in trait_cols:
                traitcolix = df.columns.get_loc(trait_col) if isinstance(trait_col, str) else trait_col
                colmask[traitcolix] = True

        # extract trait and matrix data
        trait = df.columns[colmask].to_numpy(dtype = object)
        mat = df.iloc[:,colmask].to_numpy(dtype = float)
            
        # construct output object
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

    @classmethod
    def from_csv(
            cls,
            filename: str,
            location: Union[numpy.ndarray,Real] = 0.0, 
            scale: Union[numpy.ndarray,Real] = 1.0, 
            taxa_col: Optional[Union[str,Integral]] = "taxa",
            taxa_grp_col: Optional[Union[str,Integral]] = "taxa_grp",
            trait_cols: Optional[Union[str,Sequence]] = "infer",
            sep: str = ',',
            header: int = 0,
            **kwargs: dict
        ) -> 'DenseBreedingValueMatrix':
        """
        Read a DenseBreedingValueMatrix from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.
        
        sep : str, default = ','
            CSV delimiter to use.
        
        header : int, list of int, default=0
            Row number(s) to use as the column names, and the start of the data.

        location : numpy.ndarray, Real, default = 0.0
            A ``numpy.ndarray`` of shape ``(t,)`` containing breeding value 
            locations. If given a ``Real``, create a ``numpy.ndarray`` of shape 
            ``(t,)`` filled with the provided value.

        scale : numpy.ndarray, Real, default = 1.0
            A ``numpy.ndarray`` of shape ``(t,)`` containing breeding value 
            scales. If given a ``Real``, create a ``numpy.ndarray`` of shape 
            ``(t,)`` filled with the provided value.

        taxa_col : str, Integral, None, default = "taxa"
            Name of the column from which to read taxa names.
            If of type ``str``, taxa names are read from the column named 
            defined by ``taxa_col``.
            If of type ``Integral``, taxa names are read from the column 
            number defined by ``taxa_col``.
            If ``None``, taxa names are not imported.
        
        taxa_grp_col : str, None, default = "taxa_grp"
            Name of the column from which to read taxa group names.
            If of type ``str``, taxa group names are read from the column named 
            defined by ``taxa_col``.
            If of type ``Integral``, taxa group names are read from the column 
            number defined by ``taxa_col``.
            If ``None``, taxa group names are not imported.

        trait_cols : Sequence, str, None, default = "trait"
            Names of the trait columns to which to read breeding values.
            If ``Sequence``, column names are given by the strings or integers 
            in the ``trait_cols`` Sequence.
            If ``str``, must be equal to ``"infer"``. Use remaining columns in 
            the input dataframe to load trait breeding values.
            If ``None``, do not load any trait breeding values.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : DenseBreedingValueMatrix
            A DenseBreedingValueMatrix read from a CSV file.
        """
        # read file using pandas
        df = pandas.read_csv(
            filepath_or_buffer = filename,
            sep = sep,
            header = header,
            **kwargs
        )

        # construct genetic map from pandas.DataFrame
        out = cls.from_pandas(
            df = df,
            location = location, 
            scale = scale, 
            taxa_col = taxa_col, 
            taxa_grp_col = taxa_grp_col, 
            trait_cols = trait_cols, 
            **kwargs
        )

        return out

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
            check_h5py_File_has_group(h5file, filename, groupname)    # check that group exists
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
            check_h5py_File_has_group(h5file, filename, fieldname)    # check that group exists
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
        data_dict["taxa"] = numpy.array(                        # convert taxa strings from byte to utf-8
            [s.decode("utf-8") for s in data_dict["taxa"]],
            dtype = object
        )
        data_dict["trait"] = numpy.array(                       # convert trait string from byte to utf-8
            [s.decode("utf-8") for s in data_dict["trait"]],
            dtype = object
        )
        ######################################################### create object
        gmat = cls(**data_dict)                                 # create object from read data
        return gmat

    @classmethod
    def from_numpy(
            cls, 
            mat: numpy.ndarray, 
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
        check_ndarray_ndim(mat, "mat", 2)

        # calculate location parameters
        # (n,t) -> (t,)
        location = numpy.nanmean(mat, axis = 0)

        # calculate scale parameters
        # (n,t) -> (t,)
        scale = numpy.nanstd(mat, axis = 0)

        # if scale < tolerance, set to 1.0 (do not scale)
        mask = (scale == 0.0)
        scale[mask] = 1.0

        # mean center and scale values
        # scalar / (t,) -> (t,)
        # (t,) * ( (n,t) - (t,) ) -> (n,t)
        # multiply since multiplication is faster than division for floating points
        mat = (1.0 / scale) * (mat - location) 

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



################################## Utilities ###################################
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
