"""
Module providing dense coancestry matrix implementations and associated error checking routines.
"""

__all__ = [
    "DenseCoancestryMatrix",
    "check_is_DenseCoancestryMatrix",
]

from numbers import Integral, Real
from pathlib import Path
from typing import Optional, Sequence, Union
import numpy
import warnings
import h5py
from numpy.typing import DTypeLike
import pandas

from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_pandas import check_Series_all_type, check_is_pandas_DataFrame
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group, check_h5py_File_is_readable, check_h5py_File_is_writable
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_column
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_column_index
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_column_indices
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_columns
from pybrops.core.error.error_value_pandas import check_pandas_Series_has_indices
from pybrops.core.error.error_value_pandas import check_pandas_Series_has_values
from pybrops.core.error.error_value_python import check_all_equal
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_numpy import check_ndarray_dtype
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_object
from pybrops.core.error.error_value_numpy import check_ndarray_has_values, check_ndarray_ndim
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len
from pybrops.core.error.error_type_python import check_Sequence_all_type
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_type_python import check_is_str_or_Integral
from pybrops.core.error.error_type_python import check_is_str_or_Sequence
from pybrops.core.error.error_value_python import check_str_value
from pybrops.core.mat.DenseSquareTaxaMatrix import DenseSquareTaxaMatrix
from pybrops.core.util.h5py import h5py_File_write_dict
from pybrops.popgen.cmat.CoancestryMatrix import CoancestryMatrix

class DenseCoancestryMatrix(
        DenseSquareTaxaMatrix,
        CoancestryMatrix,
    ):
    """
    A concrete class for dense coancestry matrices. Coancestry matrices are square.

    The purpose of this concrete class is to implement functionality for:
        1) Dense coancestry matrix value calculation.
        2) Dense coancestry matrix value access.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseCoancestryMatrix class.

        Parameters
        ----------
        mat : numpy.ndarray
            Coancestry matrix of shape ``(n,n)``.

            Where:

            - ``n`` is the number of taxa.
        taxa : numpy.ndarray
            Taxa names.
        taxa_grp : numpy.ndarray
            Taxa groupings.
        kwargs : dict
            Additional keyword arguments.
        """
        # call constructor for DenseTaxaMatrix
        super(DenseCoancestryMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    ############################ Object Properties #############################

    ############## Coancestry Data Properites ##############
    @DenseSquareTaxaMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        check_is_ndarray(value, "mat")
        check_all_equal(value.shape, "mat.shape")
        check_ndarray_dtype(value, "mat", numpy.float64)
        check_ndarray_ndim(value, "mat", 2)
        self._mat = value

    ################# Taxa Data Properites #################
    @DenseSquareTaxaMatrix.taxa.setter
    def taxa(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa label array"""
        if value is not None:
            check_is_ndarray(value, "taxa")
            check_ndarray_dtype_is_object(value, "taxa")
            check_ndarray_ndim(value, "taxa", 1)
            check_ndarray_axis_len(value, "taxa", 0, self._mat.shape[0])
        self._taxa = value

    @DenseSquareTaxaMatrix.taxa_grp.setter
    def taxa_grp(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group label array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp")
            check_ndarray_dtype(value, "taxa_grp", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp", 1)
            check_ndarray_axis_len(value, "taxa_grp", 0, self._mat.shape[0])
        self._taxa_grp = value

    ############### Taxa Metadata Properites ###############
    @DenseSquareTaxaMatrix.taxa_grp_name.setter
    def taxa_grp_name(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group name array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp_name")
            check_ndarray_dtype(value, "taxa_grp_name", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_name", 1)
        self._taxa_grp_name = value

    @DenseSquareTaxaMatrix.taxa_grp_stix.setter
    def taxa_grp_stix(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group start indices array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp_stix")
            check_ndarray_dtype(value, "taxa_grp_stix", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_stix", 1)
        self._taxa_grp_stix = value

    @DenseSquareTaxaMatrix.taxa_grp_spix.setter
    def taxa_grp_spix(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group stop indices array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp_spix")
            check_ndarray_dtype(value, "taxa_grp_spix", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_spix", 1)
        self._taxa_grp_spix = value

    @DenseSquareTaxaMatrix.taxa_grp_len.setter
    def taxa_grp_len(self, value: Union[numpy.ndarray,None]) -> None:
        if value is not None:
            check_is_ndarray(value, "taxa_grp_len")
            check_ndarray_dtype(value, "taxa_grp_len", numpy.int64)
            check_ndarray_ndim(value, "taxa_grp_len", 1)
        self._taxa_grp_len = value

    ############################## Object Methods ##############################

    ################## Matrix conversion ###################
    def mat_asformat(self, format: str) -> numpy.ndarray:
        """
        Get matrix in a specific format.
        
        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        
        Returns
        -------
        out : numpy.ndarray
            Matrix in the desired output format.
        """
        # check values
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")

        # process string
        if format == "coancestry":
            return self._mat.copy()
        elif format == "kinship":
            return 0.5 * self._mat

    ############## Coancestry/kinship Methods ##############
    def coancestry(
            self, 
            *args: tuple, 
            **kwargs: dict
        ) -> Real:
        """
        Retrieve the coancestry between individuals.

        Parameters
        ----------
        args : tuple
            A tuple of matrix indices to access the coancestry.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : Real
            The coancestry between individuals.
        """
        return self._mat[args]

    def kinship(
            self, 
            *args: tuple, 
            **kwargs: dict
        ) -> Real:
        """
        Retrieve the kinship between individuals.

        Parameters
        ----------
        args : tuple
            A tuple of matrix indices to access the kinship.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : Real
            The kinship between individuals.
        """
        return 0.5 * self._mat[args]

    def is_positive_semidefinite(
            self, 
            eigvaltol: Real = 2e-14
        ) -> bool:
        """
        Determine whether the coancestry matrix is positive semidefinite.
        
        Parameters
        ----------
        eigvaltol : Real
            Eigenvalue tolerance for determining positive semidefiniteness.
            If provided eigenvalue tolerance is less than zero, the tolerance
            is set to 0.0.

        Returns
        -------
        out : bool
            Whether the coancestry matrix is positive semidefinite.
        """
        # if tolerance is less than zero, set to zero
        if eigvaltol < 0.0:
            eigvaltol = 0.0
        return numpy.all(numpy.linalg.eigvals(self._mat) >= eigvaltol)
    
    def apply_jitter(
            self, 
            eigvaltol: Real = 2e-14, 
            minjitter: Real = 1e-10, 
            maxjitter: Real = 1e-6, 
            nattempt: Integral = 100
        ) -> bool:
        """
        Add a random jitter value to the diagonal of the coancestry matrix until 
        all eigenvalues exceed the provided eigenvalue tolerance.
        This ensures that a matrix can be decomposed using the Cholesky decomposition.
        This routine attempts to apply a jitter 100 times before giving up.

        Parameters
        ----------
        eigvaltol : float
            Eigenvalue tolerance for determining positive semidefiniteness.
            If provided eigenvalue tolerance is less than zero, the tolerance
            is set to 0.0.
        minjitter : float
            Minimum jitter value applied to a diagonal element.
        maxjitter : float
            Maximum jitter value applied to a diagonal element.
        nattempt : int
            Number of jitter application attempts.
        
        Returns
        -------
        out : bool
            Whether the jitter was successfully applied.
        """
        diagix = numpy.diag_indices_from(self._mat)
        mat_diag_old = self._mat[diagix].copy()
        counter = 0
        is_not_posdef = not self.is_positive_semidefinite(eigvaltol)

        # attempt to apply jitter in nattempt or less attempts
        while (is_not_posdef) and (counter < nattempt):
            self._mat[diagix] = mat_diag_old + numpy.random.uniform(minjitter, maxjitter, len(mat_diag_old))
            is_not_posdef = not self.is_positive_semidefinite(eigvaltol)
            counter += 1
        
        # if we still weren't able to find appropriate jitter, then give old diagonals and warn
        if is_not_posdef:
            self._mat[diagix] = mat_diag_old
            warnings.warn("Unable to successfully apply jitter to meet eigenvalue tolerance")
            return False
        
        return True

    def max_inbreeding(
            self,
            format: str = "coancestry"
        ) -> Real:
        """
        Calculate the maximum attainable inbreeding after one generation for 
        the coancestry matrix. For coancestry, this is equivalent to:
        
        ..math:
            \\max(\\mathrm{diag}(\\mathbf{G}))

        or for kinship, the equivalent is:

        ..math:
            \\max(\\mathrm{diag}(\\mathbf{K}))

        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        
        Returns
        -------
        out : Real
            The maximum attainable inbreeding after one generation.
        """
        # check values
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")

        # get output in coancestry format
        out = self.mat.diagonal().max()

        # if requested format is kinship, then convert coancestry to kinship
        if format == "kinship":
            out = 0.5 * out
        
        return out

    def min_inbreeding(
            self,
            format: str = "coancestry"
        ) -> Real:
        """
        Calculate the minimum attainable inbreeding after one generation for 
        the coancestry matrix. For coancestry, this is equivalent to:
        
        ..math:
            \\frac{1}{\\mathbf{1'G1}}

        or for kinship, the equivalent is:

        ..math:
            \\frac{1}{\\mathbf{1'K1}}

        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        
        Returns
        -------
        out : Real
            The minimum attainable inbreeding after one generation.
        """
        # check values
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")

        # calculate G inverse
        Ginv = numpy.linalg.inv(self.mat)

        # calculate min inbreeding
        out = 1.0 / Ginv.sum()

        # if requested format is kinship, then convert coancestry to kinship
        if format == "kinship":
            out = 0.5 * out
        
        return out

    def inverse(
            self,
            format: str = "coancestry"
        ) -> numpy.ndarray:
        """
        Calculate the inverse of the coancestry matrix.

        Parameters
        ----------
        format : str
            Desired matrix type on which to calculate the inverse. 
            Options are "coancestry", "kinship".

        Returns
        -------
        out : numpy.ndarray
            Inverse of the coancestry or kinship matrix.
        """
        # type checks
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")

        # invert matrix
        if format == "coancestry":
            out = numpy.linalg.inv(self._mat)
        elif format == "kinship":
            out = numpy.linalg.inv(0.5 * self._mat)
        
        return out

    ############## Matrix summary statistics ###############
    def max(
            self,
            format: str = "coancestry",
            axis: Union[int,tuple,None] = None
    ) -> Union[Real,numpy.ndarray]:
        """
        Calculate the maximum coancestry or kinship for the CoancestryMatrix
        along a specified axis.

        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        axis : int, tuple of ints, None
            Axis along which to find the maximum value.
        
        Returns
        -------
        out : Real, numpy.ndarray
            Maximum coancestry or kinship for the CoancestryMatrix along the 
            specified axis.
        """
        # type checks
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")
        
        # calculate mean values
        out = self._mat.max(axis = axis)

        # process string and apply transformations
        if format == "kinship":
            out *= 0.5
        
        return out

    def mean(
            self,
            format: str = "coancestry",
            axis: Union[int,tuple,None] = None,
            dtype: Optional[DTypeLike] = None
        ) -> Real:
        """
        Calculate the mean coancestry or kinship for the CoancestryMatrix.

        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        axis : int, tuple of ints, None
            Axis along which to find the mean value.
        dtype : DTypeLike, None
            Type to use in computing the mean. If ``None`` use the native 
            float type.

        Returns
        -------
        out : Real
            Mean coancestry or kinship for the CoancestryMatrix.
        """
        # type checks
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")
        
        # calculate mean values
        out = self._mat.mean(axis = axis, dtype = dtype)

        # process string and apply transformations
        if format == "kinship":
            out *= 0.5
        
        return out

    def min(
            self,
            format: str = "coancestry",
            axis: Union[int,tuple,None] = None
    ) -> Union[Real,numpy.ndarray]:
        """
        Calculate the minimum coancestry or kinship for the CoancestryMatrix
        along a specified axis.

        Parameters
        ----------
        format : str
            Desired output format. Options are "coancestry", "kinship".
        axis : int, tuple of ints, None
            Axis along which to find the minimum value.
        
        Returns
        -------
        out : Real, numpy.ndarray
            Minimum coancestry or kinship for the CoancestryMatrix along the 
            specified axis.
        """
        # type checks
        check_is_str(format, "format")
        format = format.lower()
        check_str_value(format, "format", "coancestry", "kinship")
        
        # calculate mean values
        out = self._mat.min(axis = axis)

        # process string and apply transformations
        if format == "kinship":
            out *= 0.5
        
        return out

    ###################### Matrix I/O ######################
    def to_pandas(
            self, 
            taxa_col: str = "taxa",
            taxa_grp_col: Optional[str] = "taxa_grp",
            taxa: Union[str,Sequence[Union[str,Integral]]] = "all",
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Export a ``DenseCoancestryMatrix`` to a ``pandas.DataFrame``.

        Parameters
        ----------
        taxa_col : str, default = "taxa"
            Name of the column to which to write taxa names. Cannot be ``None``.
        
        taxa_grp_col : str, None, default = "taxa_grp"
            Name of the column to which to write taxa group names.
            If ``str``, the column is given the name in ``taxa_grp_col``.
            If ``None``, the column is not exported.

        taxa : str, Sequence, default = "all"
            Name(s) of the taxa columns for which to write coancestry values.
            If ``Sequence``, export the taxa names given by the string or 
            integer value in the ``taxa`` Sequence.
            If ``str``, must be equal to ``"all"``. Export all taxa names as is.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            ``pandas.DataFrame``.
        
        Returns
        -------
        out : pandas.DataFrame
            An output dataframe.
        """
        ###
        ### type checks
        ###

        ### taxa_col
        check_is_str(taxa_col, "taxa_col")

        ### taxa_grp_col
        if taxa_grp_col is not None:
            check_is_str(taxa_grp_col, "taxa_grp_col")

        ### taxa
        if isinstance(taxa, str):
            check_str_value(taxa, "taxa", "all")
        elif isinstance(taxa, Sequence):
            # if no taxa labels present, check that taxa names are all indices 
            if self.taxa is None:
                check_Sequence_all_type(taxa, "taxa", Integral)
            # if taxa labels present, check that taxa names are all str or indices
            # for str values, make sure they are in the taxa labels
            else:
                check_Sequence_all_type(taxa, "taxa", (str,Integral))
                taxa_str = tuple(e for e in taxa if isinstance(e,str))
                check_ndarray_has_values(self.taxa, "CoancestryMatrix.taxa", *taxa_str)
        else:
            check_is_str_or_Sequence(taxa, "taxa")

        ###
        ### Vector extraction
        ###

        # construct taxa selection indices
        taxaix = None
        # if taxa is string, get all indices
        if isinstance(taxa, str) and (taxa == "all"):
            taxaix = numpy.arange(self.ntaxa)
        # else taxa is Sequence, get select indices
        else:
            # if no taxa labels, then all taxa elements are indices
            if self.taxa is None:
                taxaix = numpy.array(taxa, dtype = int)
            # else taxa labels, then convert string to indices
            else:
                # get first instance of string or index
                taxaix = numpy.array(
                    [numpy.argmax(self.taxa == t) if isinstance(t,str) else t for t in taxa], 
                    dtype = int
                )

        # extract the central matrix
        mat_data = None
        # if taxa is string, get all values
        if isinstance(taxa, str) and (taxa == "all"):
            mat_data = self.mat
        # else taxa is Sequence, get select indices
        else:
            mat_data = self.mat[taxaix,:][:,taxaix]

        # extract taxa labels
        taxa_data = None
        # if taxa is string, get all taxa names
        if isinstance(taxa, str) and (taxa == "all"):
            if self.taxa is None:
                taxa_data = [str(i) for i in range(self.ntaxa)]
            else:
                taxa_data = [str(e) for e in self.taxa]
        # else taxa is Sequence, get select taxa names
        else:
            if self.taxa is None:
                taxa_data = [str(i) for i in taxaix]
            else:
                taxa_data = [str(e) for e in self.taxa[taxaix]]

        # extract taxa group labels
        taxa_grp_data = None
        if taxa_grp_col is not None:
            # if taxa is string, get all taxa group labels
            if isinstance(taxa, str) and (taxa == "all"):
                if self.taxa_grp is None:
                    taxa_grp_data = numpy.repeat(None, self.ntaxa)
                else:
                    taxa_grp_data = self.taxa_grp
            # else taxa is Sequence, get select taxa group labels
            else:
                if self.taxa_grp is None:
                    taxa_grp_data = numpy.repeat(None, len(taxaix))
                else:
                    taxa_grp_data = self.taxa_grp[taxaix]

        ###
        ### DataFrame construction
        ###

        # make dictionary for labels
        labels_dict = {taxa_col: taxa_data}
        if taxa_grp_col is not None:
            labels_dict[taxa_grp_col] = taxa_grp_data
        
        # make dataframe for labels
        labels_df = pandas.DataFrame(labels_dict)

        # make dataframe for values
        values_df = pandas.DataFrame(mat_data, columns = taxa_data)

        # concatenate columns
        out = pandas.concat([labels_df,values_df], axis = 1)

        return out

    def to_csv(
            self, 
            filename: str,
            taxa_col: str = "taxa",
            taxa_grp_col: Optional[str] = "taxa_grp",
            taxa: Union[str,Sequence[Union[str,Integral]]] = "all",
            sep: str = ',', 
            header: bool = True, 
            index: bool = False, 
            **kwargs: dict
        ) -> None:
        """
        Write an object to a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name to which to write.

        taxa_col : str, default = "taxa"
            Name of the column to which to write taxa names. Cannot be ``None``.
        
        taxa_grp_col : str, None, default = "taxa_grp"
            Name of the column to which to write taxa group names.
            If ``str``, the column is given the name in ``taxa_grp_col``.
            If ``None``, the column is not exported.

        taxa : str, Sequence, default = "all"
            Name(s) of the taxa columns for which to write coancestry values.
            If ``Sequence``, export the taxa names given by the string or 
            integer value in the ``taxa`` Sequence.
            If ``str``, must be equal to ``"all"``. Export all taxa names as is.

        sep : str, default = ","
            Separator to use in the exported CSV file.
        
        header : bool, default = True
            Whether to save header names.
        
        index : bool, default = False
            Whether to save a row index in the exported CSV file.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a CSV.
        """
        # construct dataframe
        df = self.to_pandas(
            taxa_col = taxa_col,
            taxa_grp_col = taxa_grp_col,
            taxa = taxa,
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
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            overwrite: bool = True,
        ) -> None:
        """
        Write ``DenseCoancestryMatrix`` to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name to which to write. File is closed after writing.
            If ``h5py.File``, an opened HDF5 file to which to write. File is not closed after writing.

        groupname : str, None
            If ``str``, an HDF5 group name under which the ``DenseCoancestryMatrix`` data is stored.
            If ``None``, the ``DenseCoancestryMatrix`` is written to the base HDF5 group.

        overwrite : bool
            Whether to overwrite values in an HDF5 file if a field already exists.
        """
        # call super function
        super(DenseCoancestryMatrix, self).to_hdf5(
            filename  = filename,
            groupname = groupname,
            overwrite = overwrite,
        )

    ############################## Class Methods ###############################

    ###################### Matrix I/O ######################
    @classmethod
    def from_pandas(
            cls, 
            df: pandas.DataFrame,
            taxa_col: Union[str,Integral] = "taxa",
            taxa_grp_col: Optional[Union[str,Integral]] = "taxa_grp",
            taxa: Union[str,Sequence[Union[str,Integral]]] = "all",
            **kwargs: dict
        ) -> 'DenseCoancestryMatrix':
        """
        Read a ``DenseCoancestryMatrix`` from a ``pandas.DataFrame``.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.

        taxa_col : str, Integral, None, default = "taxa"
            Name of the column from which to read taxa names. Cannot be ``None``.
            This column is used to search for column names to extract coancestry values.
            Elements in this column are interpreted as strings.
            If of type ``str``, taxa names are read from the column named 
            defined by ``taxa_col``.
            If of type ``Integral``, taxa names are read from the column 
            number defined by ``taxa_col``.
        
        taxa_grp_col : str, None, default = "taxa_grp"
            Name of the column from which to read taxa group names.
            If of type ``str``, taxa group names are read from the column named 
            defined by ``taxa_col``.
            If of type ``Integral``, taxa group names are read from the column 
            number defined by ``taxa_col``.
            If ``None``, taxa group names are not imported.

        taxa : str, Sequence, default = "all"
            Name(s) of the taxa columns for which to read coancestry values.
            If ``Sequence``, read the taxa names given by the string or 
            integer value in the ``taxa`` Sequence.
            If ``str``, must be equal to ``"all"``. Import all taxa names as 
            defined in the ``taxa_col`` column.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            ``pandas.DataFrame``.

        Returns
        -------
        out : DenseCoancestryMatrix
            A ``DenseCoancestryMatrix`` read from a ``pandas.DataFrame``.
        """
        ###
        ### type checks
        ###

        ### df
        check_is_pandas_DataFrame(df, "df")

        # force conversion of column names to string
        df.columns = df.columns.astype(str)

        ### taxa_col
        if isinstance(taxa_col, str):
            check_pandas_DataFrame_has_column(df, "df", taxa_col)
        elif isinstance(taxa_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", taxa_col)
        else:
            check_is_str_or_Integral(taxa_col, "taxa_col")

        ### taxa_grp_col
        if taxa_grp_col is not None:
            if isinstance(taxa_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", taxa_grp_col)
                col = df[taxa_grp_col]
            elif isinstance(taxa_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", taxa_grp_col)
                col = df.iloc[:,taxa_grp_col]
                # if all entries are None, treat as though no taxa group column was provided
                if all(e is None for e in col):
                    taxa_grp_col = None
            else:
                check_is_str_or_Integral(taxa_grp_col, "taxa_grp_col")
        # extract the taxa series for future use
        taxa_col_pds = df[taxa_col] if isinstance(taxa_col,str) else df.iloc[:,taxa_col]

        ### taxa
        if isinstance(taxa, str):
            check_str_value(taxa, "taxa", "all")
            check_Series_all_type(taxa_col_pds, "df[taxa_col]", (str,Integral))
            taxa_str = tuple(e for e in taxa_col_pds if isinstance(e,str))
            taxa_int = tuple(e for e in taxa_col_pds if isinstance(e,Integral))
            check_pandas_DataFrame_has_columns(df, "df", *taxa_str)
            check_pandas_DataFrame_has_column_indices(df, "df", *taxa_int)
        elif isinstance(taxa, Sequence):
            check_Series_all_type(taxa, "taxa", (str,Integral))
            taxa_str = tuple(e for e in taxa if isinstance(e,str))
            taxa_int = tuple(e for e in taxa if isinstance(e,Integral))
            check_pandas_DataFrame_has_columns(df, "df", *taxa_str)
            check_pandas_DataFrame_has_column_indices(df, "df", *taxa_int)
            check_pandas_Series_has_values(taxa_col_pds, "df[taxa_col]", *taxa_str)
            check_pandas_Series_has_indices(taxa_col_pds, "df[taxa_col]", *taxa_int)
        else:
            check_is_str_or_Sequence(taxa, "taxa")

        ### 
        ### sanitize inputs
        ### 

        if taxa_grp_col is not None:
            # calculate column index
            colix = df.columns.get_loc(taxa_grp_col) if isinstance(taxa_grp_col,str) else taxa_grp_col
            # get mask of where values are NA
            mask = df.iloc[:,colix].isna()
            # convert None to pandas.NA
            df.iloc[mask,colix] = pandas.NA
            # if all entries are None, treat as though no taxa group column was provided
            if mask.all():
                taxa_grp_col = None

        ###
        ### data extraction
        ###

        # calculate row indices from the taxa column index
        rowix = None
        if isinstance(taxa, str) and (taxa == "all"):
            rowix = [i for i in range(len(taxa_col_pds))]
        elif isinstance(taxa, Sequence):
            taxa_col_pdi = pandas.Index(taxa_col_pds)
            rowix = [taxa_col_pdi.get_loc(e) if isinstance(e,str) else e for e in taxa]
        
        # calculate column indices
        colix = None
        if isinstance(taxa, str) and (taxa == "all"):
            colix = [df.columns.get_loc(str(e)) for e in taxa_col_pds]
        elif isinstance(taxa, Sequence):
            colix = [df.columns.get_loc(str(e)) for e in taxa]

        # extract dataframe rows and columns for matrix construction
        mat_data = df.iloc[rowix,:].iloc[:,colix].to_numpy(dtype = float, na_value = numpy.nan)
        taxa_data = taxa_col_pds[rowix].to_numpy(dtype = object)
        taxa_grp_data = None
        if taxa_grp_col is not None:
            taxa_grp_col_pds = df[taxa_grp_col] if isinstance(taxa_grp_col,str) else df.iloc[:,taxa_grp_col]
            taxa_grp_data = taxa_grp_col_pds[rowix].to_numpy(dtype = int, na_value = -1)

        ###
        ### construct object
        ###

        # construct object
        out = cls(
            mat = mat_data,
            taxa = taxa_data,
            taxa_grp = taxa_grp_data,
            **kwargs
        )

        return out

    @classmethod
    def from_csv(
            cls, 
            filename: str,
            sep: str = ',',
            header: int = 0,
            taxa_col: Union[str,Integral] = "taxa",
            taxa_grp_col: Optional[Union[str,Integral]] = "taxa_grp",
            taxa: Union[str,Sequence[Union[str,Integral]]] = "all",
            **kwargs: dict
        ) -> 'DenseCoancestryMatrix':
        """
        Read a ``DenseCoancestryMatrix`` from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.
        
        sep : str, default = ','
            CSV delimiter to use.
        
        header : int, list of int, default=0
            Row number(s) to use as the column names, and the start of the data.

        taxa_col : str, Integral, None, default = "taxa"
            Name of the column from which to read taxa names. Cannot be ``None``.
            This column is used to search for column names to extract coancestry values.
            Elements in this column are interpreted as strings.
            If of type ``str``, taxa names are read from the column named 
            defined by ``taxa_col``.
            If of type ``Integral``, taxa names are read from the column 
            number defined by ``taxa_col``.
        
        taxa_grp_col : str, None, default = "taxa_grp"
            Name of the column from which to read taxa group names.
            If of type ``str``, taxa group names are read from the column named 
            defined by ``taxa_col``.
            If of type ``Integral``, taxa group names are read from the column 
            number defined by ``taxa_col``.
            If ``None``, taxa group names are not imported.

        taxa : str, Sequence, default = "all"
            Name(s) of the taxa columns for which to read coancestry values.
            If ``Sequence``, read the taxa names given by the string or 
            integer value in the ``taxa`` Sequence.
            If ``str``, must be equal to ``"all"``. Import all taxa names as 
            defined in the ``taxa_col`` column.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : DenseCoancestryMatrix
            A ``DenseCoancestryMatrix`` read from a CSV file.
        """
        # type checks
        check_is_str(sep, "sep")

        # read file using pandas
        df = pandas.read_csv(
            filename,
            sep = sep,
            header = header,
            **kwargs
        )

        # construct genetic map from pandas.DataFrame
        out = cls.from_pandas(
            df = df,
            taxa_col = taxa_col, 
            taxa_grp_col = taxa_grp_col, 
            taxa = taxa, 
            **kwargs
        )

        return out

    @classmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None
        ) -> 'DenseCoancestryMatrix':
        """
        Read a ``DenseCoancestryMatrix`` from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which the ``DenseCoancestryMatrix`` data is stored.
            If ``None``, the ``DenseCoancestryMatrix`` is read from base HDF5 group.

        Returns
        -------
        out : DenseCoancestryMatrix
            A ``DenseCoancestryMatrix`` read from file.
        """
        # call super function
        return super(DenseCoancestryMatrix, cls).from_hdf5(
            filename  = filename,
            groupname = groupname,
        )



################################## Utilities ###################################
def check_is_DenseCoancestryMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseCoancestryMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseCoancestryMatrix):
        raise TypeError("variable '{0}' must be a DenseCoancestryMatrix".format(vname))
