"""
Module implementing classes and associated error checking routines for matrices
storing dense additive genic covariance estimates calculated using three-way DH
formulae.
"""

import math
from numbers import Integral
from pathlib import Path
from typing import Optional, Union
import numpy
import pandas
import h5py
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_type_python import check_is_Integral, check_is_str, check_is_str_or_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_column, check_pandas_DataFrame_has_column_index
from pybrops.core.util.arrayix import flattenix
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel, check_is_AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.model.pcvmat.DenseAdditiveProgenyGenicCovarianceMatrix import DenseAdditiveProgenyGenicCovarianceMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix

class DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix(DenseAdditiveProgenyGenicCovarianceMatrix):
    """
    A concrete class for dense additive genic covariance matrices calculated
    for two-way DH progenies.

    The purpose of this concrete class is to implement functionality for:
        1) Genic covariance estimation for two-way DH progenies.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ):
        """
        Constructor for the concrete class DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            Array used to construct the object.
        taxa : numpy.ndarray
            Taxa names.
        taxa_grp : numpy.ndarray
            Taxa groupings.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    ############################ Object Properties #############################

    ##################### Matrix Data ######################
    @DenseAdditiveProgenyGenicCovarianceMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set pointer to raw numpy.ndarray object."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim(value, "mat", 5) # (ntaxa,ntaxa,ntaxa,ntrait,ntrait)
        self._mat = value

    ############## Square Metadata Properties ##############
    @DenseAdditiveProgenyGenicCovarianceMatrix.square_taxa_axes.getter
    def square_taxa_axes(self) -> tuple:
        """Axis indices for taxa axes that are square."""
        return (0,1,2) # (recurrent, female, male)

    @DenseAdditiveProgenyGenicCovarianceMatrix.square_trait_axes.getter
    def square_trait_axes(self) -> tuple:
        return (3,4) # (trait1, trait2) covariance matrix

    #################### Trait metadata ####################

    ######## Expected parental genome contributions ########
    @DenseAdditiveProgenyGenicCovarianceMatrix.epgc.getter
    def epgc(self) -> tuple:
        """Get a tuple of the expected parental genome contributions."""
        return (0.5, 0.25, 0.25)

    ################# Parental dimensions ##################
    @property
    def nrecurrent(self) -> Integral:
        """Number of recurrent parents."""
        return self._mat.shape[self.recurrent_axis]
    
    @property
    def recurrent_axis(self) -> Integral:
        """Axis along which recurrent parents are stored."""
        return 0

    @property
    def nfemale(self) -> Integral:
        """Number of female parents."""
        return self._mat.shape[self.female_axis]
    
    @property
    def female_axis(self) -> Integral:
        """Axis along which female parents are stored."""
        return 1
    
    @property
    def nmale(self) -> Integral:
        """Number of male parents."""
        return self._mat.shape[self.male_axis]
    
    @property
    def male_axis(self) -> Integral:
        """Axis along which male parents are stored."""
        return 2

    ############################## Object Methods ##############################

    ###################### Matrix I/O ######################
    # TODO: make exporting of specific female, male, trait values, not just all
    def to_pandas(
            self, 
            recurrent_col: str = "recurrent",
            recurrent_grp_col: Optional[str] = "recurrent_grp",
            female_col: str = "female",
            female_grp_col: Optional[str] = "female_grp",
            male_col: str = "male",
            male_grp_col: Optional[str] = "male_grp",
            trait1_col: str = "trait1",
            trait2_col: str = "trait2",
            covariance_col: str = "covariance",
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Export a DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix to a pandas.DataFrame.

        Parameters
        ----------
        female_col : str, default = "female"
            Name of the column to which to write female taxa names.

        female_grp_col : str, None, default = "female_grp"
            Name of the column to which to write female taxa groups.

        male_col : str, default = "male"
            Name of the column to which to write male taxa names.

        male_grp_col : str, None, default = "male_grp"
            Name of the column to which to write male taxa groups.

        recurrent_col : str, default = "recurrent"
            Name of the column to which to write recurrent taxa names.

        recurrent_grp_col : str, None, default = "recurrent_grp"
            Name of the column to which to write recurrent taxa groups.

        trait1_col : str, default = "trait"
            Name of the column to which to write trait taxa names.

        covariance_col : str, default = "covariance"
            Name of the column to which to write covariance taxa names.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            pandas.DataFrame.
        
        Returns
        -------
        out : pandas.DataFrame
            An output dataframe.
        """
        # type checks
        check_is_str(female_col, "female_col")
        if female_grp_col is not None:
            check_is_str(female_grp_col, "female_grp_col")

        check_is_str(male_col, "male_col")
        if male_grp_col is not None:
            check_is_str(male_grp_col, "male_grp_col")
            
        check_is_str(recurrent_col, "recurrent_col")
        if recurrent_grp_col is not None:
            check_is_str(recurrent_grp_col, "recurrent_grp_col")
            
        check_is_str(trait1_col, "trait1_col")
        check_is_str(trait2_col, "trait2_col")
        check_is_str(covariance_col, "covariance_col")

        # get names for taxa
        if self.taxa is None:
            nzero = math.ceil(math.log10(self.ntaxa))+1
            taxa = numpy.array(["Taxon"+str(e).zfill(nzero) for e in range(self.ntaxa)], dtype = object)
        else:
            taxa = self.taxa

        # get names for traits
        if self.trait is None:
            nzero = math.ceil(math.log10(self.ntrait))+1
            trait = numpy.array(["Trait"+str(e).zfill(nzero) for e in range(self.ntrait)], dtype = object)
        else:
            trait = self.trait

        # calculate flattened array and corresponding axis indices
        flatmat, (recurix, femaleix, maleix, trait1ix, trait2ix) = flattenix(self.mat)

        # make dictionary to store output columns in specific column ordering
        out_dict = {}
        # add recurrent parent values
        out_dict.update({recurrent_col: taxa[recurix]})
        if recurrent_grp_col is not None:
            values = None if self.taxa_grp is None else self.taxa_grp[recurix]
            out_dict.update({recurrent_grp_col: values})
        # add female parent values
        out_dict.update({female_col: taxa[femaleix]})
        if female_grp_col is not None:
            values = None if self.taxa_grp is None else self.taxa_grp[femaleix]
            out_dict.update({female_grp_col: values})
        # add male parent values
        out_dict.update({male_col: taxa[maleix]})
        if male_grp_col is not None:
            values = None if self.taxa_grp is None else self.taxa_grp[maleix]
            out_dict.update({male_grp_col: values})
        # add trait and covariance values
        out_dict.update({trait1_col: trait[trait1ix]})
        out_dict.update({trait2_col: trait[trait2ix]})
        out_dict.update({covariance_col: flatmat})

        # create a pandas DataFrame from the data
        out = pandas.DataFrame(out_dict, **kwargs)

        return out

    def to_csv(
            self, 
            filename: str,
            recurrent_col: str = "recurrent",
            recurrent_grp_col: Optional[str] = "recurrent_grp",
            female_col: str = "female",
            female_grp_col: Optional[str] = "female_grp",
            male_col: str = "male",
            male_grp_col: Optional[str] = "male_grp",
            trait1_col: str = "trait1",
            trait2_col: str = "trait2",
            covariance_col: str = "covariance",
            sep: str = ',', 
            header: bool = True, 
            index: bool = False, 
            **kwargs: dict
        ) -> None:
        """
        Write a DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix to a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name to which to write.
        
        female_col : str, default = "female"
            Name of the column to which to write female taxa names.

        female_grp_col : str, None, default = "female_grp"
            Name of the column to which to write female taxa groups.

        male_col : str, default = "male"
            Name of the column to which to write male taxa names.

        male_grp_col : str, None, default = "male_grp"
            Name of the column to which to write male taxa groups.

        trait1_col : str, default = "trait"
            Name of the column to which to write trait taxa names.

        covariance_col : str, default = "covariance"
            Name of the column to which to write covariance taxa names.

        sep : str, default = ","
            Separator to use in the exported CSV file.
        
        header : bool, default = True
            Whether to save header names.
        
        index : bool, default = False
            Whether to save a row index in the exported CSV file.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a CSV.
        """
        # convert DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix to pandas.DataFrame
        df = self.to_pandas(
            recurrent_col = recurrent_col,
            recurrent_grp_col = recurrent_grp_col,
            female_col = female_col,
            female_grp_col = female_grp_col,
            male_col = male_col,
            male_grp_col = male_grp_col,
            trait1_col = trait1_col,
            trait2_col = trait2_col,
            covariance_col = covariance_col,
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
        Write ``DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix`` to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            HDF5 file name or HDF5 file stream to which to write.

        groupname : str, None
            If ``str``, an HDF5 group name under which the ``DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix`` data is stored.
            If ``None``, the ``DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix`` is written to the base HDF5 group.
        
        overwrite : bool
            Whether to overwrite data fields if they are present in the HDF5 file.
        """
        # call super function
        super(DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix, self).to_hdf5(
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
            recurrent_col: Union[str,Integral] = "recurrent",
            recurrent_grp_col: Optional[Union[str,Integral]] = "recurrent_grp",
            female_col: Union[str,Integral] = "female",
            female_grp_col: Optional[Union[str,Integral]] = "female_grp",
            male_col: Union[str,Integral] = "male",
            male_grp_col: Optional[Union[str,Integral]] = "male_grp",
            trait1_col: Union[str,Integral] = "trait1",
            trait2_col: Union[str,Integral] = "trait2",
            covariance_col: Union[str,Integral] = "covariance",
            **kwargs: dict
        ) -> 'DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix':
        """
        Read a ``DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix`` from a ``pandas.DataFrame``.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.

        recurrent_col : str, Integral, default = "recurrent"
            Name or index of the column from which to read recurrent taxa names.

        recurrent_grp_col : str, Integral, None, default = "recurrent"
            Name or index of the column from which to read recurrent taxa groups.

        female_col : str, Integral, default = "female"
            Name or index of the column from which to read female taxa names.

        female_grp_col : str, Integral, None, default = "female"
            Name or index of the column from which to read female taxa groups.

        male_col : str, Integral, None, default = "male"
            Name or index of the column from which to read male taxa names.

        male_grp_col : str, Integral, None, default = "male"
            Name or index of the column from which to read male taxa names.

        trait1_col : str, Integral, default = "trait"
            Name or index of the column from which to read trait taxa names.

        covariance_col : str, Integral, default = "covariance"
            Name or index of the column from which to read covariance taxa names.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            ``pandas.DataFrame``.
        
        Returns
        -------
        out : DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix
            A ``DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix`` read from a ``pandas.DataFrame``.
        """
        ### type checks and get df indices

        # df
        check_is_pandas_DataFrame(df, "df")

        # recurrent_col
        if isinstance(recurrent_col, str):
            check_pandas_DataFrame_has_column(df, "df", recurrent_col)
            recurrent_colix = df.columns.get_loc(recurrent_col)
        elif isinstance(recurrent_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", recurrent_col)
            recurrent_colix = recurrent_col
        else:
            check_is_str_or_Integral(recurrent_col, "recurrent_col")
        
        # recurrent_grp_col
        if recurrent_grp_col is not None:
            if isinstance(recurrent_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", recurrent_grp_col)
                recurrent_grp_colix = df.columns.get_loc(recurrent_grp_col)
            elif isinstance(recurrent_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", recurrent_grp_col)
                recurrent_grp_colix = recurrent_grp_col
            else:
                check_is_str_or_Integral(recurrent_grp_col, "recurrent_grp_col")
        
        # female_col
        if isinstance(female_col, str):
            check_pandas_DataFrame_has_column(df, "df", female_col)
            female_colix = df.columns.get_loc(female_col)
        elif isinstance(female_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", female_col)
            female_colix = female_col
        else:
            check_is_str_or_Integral(female_col, "female_col")
        
        # female_grp_col
        if female_grp_col is not None:
            if isinstance(female_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", female_grp_col)
                female_grp_colix = df.columns.get_loc(female_grp_col)
            elif isinstance(female_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", female_grp_col)
                female_grp_colix = female_grp_col
            else:
                check_is_str_or_Integral(female_grp_col, "female_grp_col")
        
        # male_col
        if isinstance(male_col, str):
            check_pandas_DataFrame_has_column(df, "df", male_col)
            male_colix = df.columns.get_loc(male_col)
        elif isinstance(male_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", male_col)
            male_colix = male_col
        else:
            check_is_str_or_Integral(male_col, "male_col")
        
        # male_grp_col
        if male_grp_col is not None:
            if isinstance(male_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", male_grp_col)
                male_grp_colix = df.columns.get_loc(male_grp_col)
            elif isinstance(male_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", male_grp_col)
                male_grp_colix = male_grp_col
            else:
                check_is_str_or_Integral(male_grp_col, "male_grp_col")
        
        # trait1_col
        if isinstance(trait1_col, str):
            check_pandas_DataFrame_has_column(df, "df", trait1_col)
            trait1_colix = df.columns.get_loc(trait1_col)
        elif isinstance(trait1_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", trait1_col)
            trait1_colix = trait1_col
        else:
            check_is_str_or_Integral(trait1_col, "trait1_col")
        
        # trait2_col
        if isinstance(trait2_col, str):
            check_pandas_DataFrame_has_column(df, "df", trait2_col)
            trait2_colix = df.columns.get_loc(trait2_col)
        elif isinstance(trait2_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", trait2_col)
            trait2_colix = trait2_col
        else:
            check_is_str_or_Integral(trait2_col, "trait2_col")
        
        # covariance_col
        if isinstance(covariance_col, str):
            check_pandas_DataFrame_has_column(df, "df", covariance_col)
            covariance_colix = df.columns.get_loc(covariance_col)
        elif isinstance(covariance_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", covariance_col)
            covariance_colix = covariance_col
        else:
            check_is_str_or_Integral(covariance_col, "covariance_col")
        
        # get required data columns (type numpy.ndarray)
        recurrent_data   = df.iloc[:,recurrent_colix ].to_numpy(dtype = object)
        female_data      = df.iloc[:,female_colix    ].to_numpy(dtype = object)
        male_data        = df.iloc[:,male_colix      ].to_numpy(dtype = object)
        trait1_data      = df.iloc[:,trait1_colix    ].to_numpy(dtype = object)
        trait2_data      = df.iloc[:,trait2_colix    ].to_numpy(dtype = object)
        covariance_data  = df.iloc[:,covariance_colix].to_numpy(dtype = float)

        # get unique recurrent, female, male taxa (type numpy.ndarray)
        recurrent_taxa, recurrent_taxaix = numpy.unique(recurrent_data, return_index = True)
        female_taxa, female_taxaix = numpy.unique(female_data, return_index = True)
        male_taxa, male_taxaix     = numpy.unique(male_data, return_index = True)

        # combine recurrent, female, male taxa names (type numpy.ndarray)
        taxa = numpy.union1d(recurrent_taxa, numpy.union1d(female_taxa, male_taxa))

        # allocate female, male indices
        recurix  = numpy.empty(len(recurrent_data), dtype = int)
        femaleix = numpy.empty(len(female_data), dtype = int)
        maleix   = numpy.empty(len(male_data  ), dtype = int)

        # calculate female, male indices
        for i,taxon in enumerate(taxa):
            recurix[recurrent_data == taxon] = i
            femaleix[female_data == taxon] = i
            maleix[male_data == taxon] = i
        
        # calculate unique trait values, trait indices
        trait1, trait1ix = numpy.unique(trait1_data, return_inverse = True)
        trait2, trait2ix = numpy.unique(trait2_data, return_inverse = True)

        trait = numpy.union1d(trait1, trait2)

        # get optional taxa group data
        taxa_grp = None
        if recurrent_grp_col is not None:
            if taxa_grp is None:
                taxa_grp = numpy.empty(len(taxa), dtype = int)
            recurrent_grp_data = df.iloc[:,recurrent_grp_colix].to_numpy(dtype = int)
            for i,taxon in enumerate(taxa):
                if taxon in recurrent_taxa:
                    ix = recurrent_taxaix[recurrent_taxa == taxon][0]
                    taxa_grp[i] = recurrent_grp_data[ix]
        if female_grp_col is not None:
            if taxa_grp is None:
                taxa_grp = numpy.empty(len(taxa), dtype = int)
            female_grp_data = df.iloc[:,female_grp_colix].to_numpy(dtype = int)
            for i,taxon in enumerate(taxa):
                if taxon in female_taxa:
                    ix = female_taxaix[female_taxa == taxon][0]
                    taxa_grp[i] = female_grp_data[ix]
        if male_grp_col is not None:
            if taxa_grp is None:
                taxa_grp = numpy.empty(len(taxa), dtype = int)
            male_grp_data = df.iloc[:,male_grp_colix].to_numpy(dtype = int)
            for i,taxon in enumerate(taxa):
                if taxon in male_taxa:
                    ix = male_taxaix[male_taxa == taxon][0]
                    taxa_grp[i] = male_grp_data[ix]

        # get array dimensions
        nrecurrent = len(taxa)
        nfemale = len(taxa)
        nmale = len(taxa)
        ntrait = len(trait)

        # allocate NaN array for covariance matrix
        mat = numpy.full((nrecurrent,nfemale,nmale,ntrait,ntrait), numpy.nan, dtype = float)

        # overwrite NaN values with covariance values
        mat[recurix,femaleix,maleix,trait1ix,trait2ix] = covariance_data

        # construct an object
        out = cls(
            mat = mat, 
            taxa = taxa, 
            taxa_grp = taxa_grp, 
            trait = trait, 
        )

        return out

    @classmethod
    def from_csv(
            cls, 
            filename: str, 
            recurrent_col: Union[str,Integral] = "recurrent",
            recurrent_grp_col: Optional[Union[str,Integral]] = "recurrent_grp",
            female_col: Union[str,Integral] = "female",
            female_grp_col: Optional[Union[str,Integral]] = "female_grp",
            male_col: Union[str,Integral] = "male",
            male_grp_col: Optional[Union[str,Integral]] = "male_grp",
            trait1_col: Union[str,Integral] = "trait1",
            trait2_col: Union[str,Integral] = "trait2",
            covariance_col: Union[str,Integral] = "covariance",
            sep: str = ',',
            header: int = 0,
            **kwargs: dict
        ) -> 'DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix':
        """
        Read an object from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.

        recurrent_col : str, Integral, default = "recurrent"
            Name or index of the column from which to read recurrent taxa names.

        recurrent_grp_col : str, Integral, None, default = "recurrent"
            Name or index of the column from which to read recurrent taxa groups.

        female_col : str, Integral, default = "female"
            Name or index of the column from which to read female taxa names.

        female_grp_col : str, Integral, None, default = "female"
            Name or index of the column from which to read female taxa groups.

        male_col : str, Integral, None, default = "male"
            Name or index of the column from which to read male taxa names.

        male_grp_col : str, Integral, None, default = "male"
            Name or index of the column from which to read male taxa names.

        trait1_col : str, Integral, default = "trait"
            Name or index of the column from which to read trait taxa names.

        covariance_col : str, Integral, default = "covariance"
            Name or index of the column from which to read covariance taxa names.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : CSVInputOutput
            An object read from a CSV file.
        """
        df = pandas.read_csv(
            filepath_or_buffer = filename,
            sep = sep,
            header = header,
            **kwargs
        )

        out = cls.from_pandas(
            df = df,
            recurrent_col = recurrent_col,
            recurrent_grp_col = recurrent_grp_col,
            female_col = female_col,
            female_grp_col = female_grp_col,
            male_col = male_col,
            male_grp_col = male_grp_col,
            trait1_col = trait1_col,
            trait2_col = trait2_col,
            covariance_col = covariance_col,
        )

        return out
    
    @classmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None
        ) -> 'DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix':
        """
        Read ``DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix`` from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which ``DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix`` data is stored.
            If ``None``, ``DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix`` is read from base HDF5 group.

        Returns
        -------
        out : DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix
            A ``DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix`` read from file.
        """
        # call super function
        return super(DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix, cls).from_hdf5(
            filename  = filename,
            groupname = groupname,
        )

    ############# Matrix Factory Class Methods #############
    # TODO: provide support for non-linear models
    @classmethod
    def from_gmod(
            cls, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nprogeny: int,
            **kwargs: dict
        ) -> 'DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix':
        """
        Estimate genic variances from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genic variances.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genic
            covariance.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : ProgenyGenicCovarianceMatrix
            A matrix of genic covariance estimations.
        """
        # type checks
        check_is_GenomicModel(gmod, "gmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_Integral(nprogeny, "nprogeny")

        # if genomic model is an additive linear genomic model, then use specialized routine
        if isinstance(gmod, AdditiveLinearGenomicModel):
            return cls.from_algmod(
                algmod = gmod, 
                pgmat = pgmat, 
                nprogeny = nprogeny, 
                **kwargs
            )
        # otherwise raise error since non-linear support hasn't been implemented yet
        else:
            raise NotImplementedError("support for non-linear models not implemented yet")

    @classmethod
    def from_algmod(
            cls, 
            algmod: AdditiveLinearGenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nprogeny: int, 
            mem: int
        ) -> 'DenseThreeWayDHAdditiveProgenyGenicCovarianceMatrix':
        """
        Estimate genic variances from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genic variances.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genic
            covariance.
        mem : int
            Memory chunk size to use during matrix operations.

        Returns
        -------
        out : ProgenyGenicCovarianceMatrix
            A matrix of additive genic covariance estimations.
        """
        # type checks
        check_is_AdditiveLinearGenomicModel(algmod, "algmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")

        # gather shapes of data input
        ntrait = algmod.ntrait                  # number of traits (t)
        ntaxa = pgmat.ntaxa                     # number of individuals (n)
        ploidy = pgmat.ploidy                   # ploidy level (scalar)
        epgc = (0.5,0.25,0.25)                  # parental contributions

        # gather allele frequencies within each taxon
        tafreq = pgmat.tafreq()                 # (n,p) allele frequencies within taxon
        u = algmod.u_a                          # (p,t) marker effect coefficients
        
        # allocate a square matrix for each pairwise covariance
        cov_a = numpy.empty(
            (ntaxa,ntaxa,ntrait,ntrait),               # (n,n,t) covariance matrix
            dtype = float
        )

        # for each mate pair (including selfs)
        for recurr in range(0,ntaxa):
            for female in range(0,ntaxa):
                for male in range(0,female):
                    # calculate the cross allele frequency
                    # (n,p)[(3,),:] -> (3,p)
                    # (3,) . (3,p) -> (p,)
                    p = numpy.dot(epgc, tafreq[(recurr,female,male),:])

                    # calculate the diagonal matrix elements as p(1-p)
                    # scalar - (p,1) -> (p,1)
                    # (p,1) * (p,1) -> (p,t)
                    d = p[:,None] * (1.0 - p[:,None])

                    # calculate the covariance
                    # (p,t) * (p,1) -> (p,t)
                    # (p,t).T -> (t,p)
                    # (t,p) @ (p,t) -> (t,t)
                    cov = ((u * d).T) @ u

                    # multiply by the ploidy
                    cov *= ploidy

                    # store in matrix and copy to lower since matrix is symmetrical
                    cov_a[recurr,female,male,:,:] = cov
                    cov_a[recurr,male,female,:,:] = cov

        # construct output
        out = cls(
            mat = cov_a,
            taxa = pgmat.taxa,
            taxa_grp = pgmat.taxa_grp,
            trait = algmod.trait
        )

        return out
