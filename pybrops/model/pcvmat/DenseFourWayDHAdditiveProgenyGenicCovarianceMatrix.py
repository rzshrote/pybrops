"""
Module implementing classes and associated error checking routines for matrices
storing dense additive genic covariance estimates calculated using four-way DH
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

class DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix(DenseAdditiveProgenyGenicCovarianceMatrix):
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
        Constructor for the concrete class DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix.

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
        super(DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    ########### Miscellaneous special functions ############
    def __repr__(
            self
        ) -> str:
        """
        Return repr(self).
        
        Returns
        -------
        out : str
            A representation of the object.
        """
        return "<{0} of shape (nfemale2 = {1}, nmale2 = {2}, nfemale1 = {3}, nmale1 = {4}, ntrait = {5}, ntrait = {6}) at {7}>".format(
            type(self).__name__,
            self.nfemale2,
            self.nmale2,
            self.nfemale1,
            self.nmale1,
            self.ntrait,
            self.ntrait,
            hex(id(self)),
        )

    ############################ Object Properties #############################

    ##################### Matrix Data ######################
    @DenseAdditiveProgenyGenicCovarianceMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set pointer to raw numpy.ndarray object."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim(value, "mat", 6) # (ntaxa,ntaxa,ntaxa,ntaxa,ntrait,ntrait)
        self._mat = value

    ############## Square Metadata Properties ##############
    @DenseAdditiveProgenyGenicCovarianceMatrix.square_taxa_axes.getter
    def square_taxa_axes(self) -> tuple:
        """Axis indices for taxa axes that are square."""
        return (0,1,2,3) # (female2, male2, female1, male1)

    @DenseAdditiveProgenyGenicCovarianceMatrix.square_trait_axes.getter
    def square_trait_axes(self) -> tuple:
        return (4,5) # (trait1, trait2) covariance matrix

    #################### Trait metadata ####################

    ######## Expected parental genome contributions ########
    @DenseAdditiveProgenyGenicCovarianceMatrix.epgc.getter
    def epgc(self) -> tuple:
        """Get a tuple of the expected parental genome contributions."""
        return (0.25, 0.25, 0.25, 0.25)

    ################# Parental dimensions ##################
    @property
    def nfemale2(self) -> Integral:
        """Number of female 2 parents."""
        return self._mat.shape[self.female2_axis]
    
    @property
    def female2_axis(self) -> Integral:
        """Axis along which female 2 parents are stored."""
        return 0
    
    @property
    def nmale2(self) -> Integral:
        """Number of male 2 parents."""
        return self._mat.shape[self.male2_axis]
    
    @property
    def male2_axis(self) -> Integral:
        """Axis along which male 2 parents are stored."""
        return 1

    @property
    def nfemale1(self) -> Integral:
        """Number of female 1 parents."""
        return self._mat.shape[self.female1_axis]
    
    @property
    def female1_axis(self) -> Integral:
        """Axis along which female 1 parents are stored."""
        return 2
    
    @property
    def nmale1(self) -> Integral:
        """Number of male 1 parents."""
        return self._mat.shape[self.male1_axis]
    
    @property
    def male1_axis(self) -> Integral:
        """Axis along which male 1 parents are stored."""
        return 3

    ############################## Object Methods ##############################

    ###################### Matrix I/O ######################
    # TODO: make exporting of specific female, male, trait values, not just all
    def to_pandas(
            self, 
            female2_col: str = "female2",
            female2_grp_col: Optional[str] = "female2_grp",
            male2_col: str = "male2",
            male2_grp_col: Optional[str] = "male2_grp",
            female1_col: str = "female1",
            female1_grp_col: Optional[str] = "female1_grp",
            male1_col: str = "male1",
            male1_grp_col: Optional[str] = "male1_grp",
            trait1_col: str = "trait1",
            trait2_col: str = "trait2",
            covariance_col: str = "covariance",
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Export a DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix to a pandas.DataFrame.

        Parameters
        ----------
        female2_col : str, default = "female2"
            Name of the column to which to write female 2 taxa names.

        female2_grp_col : str, None, default = "female2_grp"
            Name of the column to which to write female 2 taxa groups.

        male2_col : str, default = "male2"
            Name of the column to which to write male 2 taxa names.

        male2_grp_col : str, None, default = "male2_grp"
            Name of the column to which to write male 2 taxa groups.

        female1_col : str, default = "female1"
            Name of the column to which to write female 1 taxa names.

        female1_grp_col : str, None, default = "female1_grp"
            Name of the column to which to write female 1 taxa groups.

        male1_col : str, default = "male1"
            Name of the column to which to write male 1 taxa names.

        male1_grp_col : str, None, default = "male1_grp"
            Name of the column to which to write male 1 taxa groups.

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
        check_is_str(female2_col, "female2_col")
        if female2_grp_col is not None:
            check_is_str(female2_grp_col, "female2_grp_col")
            
        check_is_str(male2_col, "male2_col")
        if male2_grp_col is not None:
            check_is_str(male2_grp_col, "male2_grp_col")
            
        check_is_str(female1_col, "female1_col")
        if female1_grp_col is not None:
            check_is_str(female1_grp_col, "female1_grp_col")

        check_is_str(male1_col, "male1_col")
        if male1_grp_col is not None:
            check_is_str(male1_grp_col, "male1_grp_col")
            
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
        flatmat, (female2ix, male2ix, female1ix, male1ix, trait1ix, trait2ix) = flattenix(self.mat)

        # make dictionary to store output columns in specific column ordering
        out_dict = {}
        # add female 2 parent values
        out_dict.update({female2_col: taxa[female2ix]})
        if female2_grp_col is not None:
            values = None if self.taxa_grp is None else self.taxa_grp[female2ix]
            out_dict.update({female2_grp_col: values})
        # add male 2 parent values
        out_dict.update({male2_col: taxa[male2ix]})
        if male2_grp_col is not None:
            values = None if self.taxa_grp is None else self.taxa_grp[male2ix]
            out_dict.update({male2_grp_col: values})
        # add female 1 parent values
        out_dict.update({female1_col: taxa[female1ix]})
        if female1_grp_col is not None:
            values = None if self.taxa_grp is None else self.taxa_grp[female1ix]
            out_dict.update({female1_grp_col: values})
        # add male 1 parent values
        out_dict.update({male1_col: taxa[male1ix]})
        if male1_grp_col is not None:
            values = None if self.taxa_grp is None else self.taxa_grp[male1ix]
            out_dict.update({male1_grp_col: values})
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
            female2_col: str = "female2",
            female2_grp_col: Optional[str] = "female2_grp",
            male2_col: str = "male2",
            male2_grp_col: Optional[str] = "male2_grp",
            female1_col: str = "female1",
            female1_grp_col: Optional[str] = "female1_grp",
            male1_col: str = "male1",
            male1_grp_col: Optional[str] = "male1_grp",
            trait1_col: str = "trait1",
            trait2_col: str = "trait2",
            covariance_col: str = "covariance",
            sep: str = ',', 
            header: bool = True, 
            index: bool = False, 
            **kwargs: dict
        ) -> None:
        """
        Write a DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix to a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name to which to write.
        
        female2_col : str, default = "female2"
            Name of the column to which to write female 2 taxa names.

        female2_grp_col : str, None, default = "female2_grp"
            Name of the column to which to write female 2 taxa groups.

        male2_col : str, default = "male2"
            Name of the column to which to write male 2 taxa names.

        male2_grp_col : str, None, default = "male2_grp"
            Name of the column to which to write male 2 taxa groups.

        female1_col : str, default = "female1"
            Name of the column to which to write female 1 taxa names.

        female1_grp_col : str, None, default = "female1_grp"
            Name of the column to which to write female 1 taxa groups.

        male1_col : str, default = "male1"
            Name of the column to which to write male 1 taxa names.

        male1_grp_col : str, None, default = "male1_grp"
            Name of the column to which to write male 1 taxa groups.

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
        # convert DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix to pandas.DataFrame
        df = self.to_pandas(
            female2_col     = female2_col,
            female2_grp_col = female2_grp_col,
            male2_col       = male2_col,
            male2_grp_col   = male2_grp_col,
            female1_col     = female1_col,
            female1_grp_col = female1_grp_col,
            male1_col       = male1_col,
            male1_grp_col   = male1_grp_col,
            trait1_col       = trait1_col,
            trait2_col       = trait2_col,
            covariance_col    = covariance_col,
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
        Write ``DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix`` to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            HDF5 file name or HDF5 file stream to which to write.

        groupname : str, None
            If ``str``, an HDF5 group name under which the ``DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix`` data is stored.
            If ``None``, the ``DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix`` is written to the base HDF5 group.
        
        overwrite : bool
            Whether to overwrite data fields if they are present in the HDF5 file.
        """
        # call super function
        super(DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix, self).to_hdf5(
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
            female2_col: Union[str,Integral] = "female2",
            female2_grp_col: Optional[Union[str,Integral]] = "female2_grp",
            male2_col: Union[str,Integral] = "male2",
            male2_grp_col: Optional[Union[str,Integral]] = "male2_grp",
            female1_col: Union[str,Integral] = "female1",
            female1_grp_col: Optional[Union[str,Integral]] = "female1_grp",
            male1_col: Union[str,Integral] = "male1",
            male1_grp_col: Optional[Union[str,Integral]] = "male1_grp",
            trait1_col: Union[str,Integral] = "trait1",
            trait2_col: Union[str,Integral] = "trait2",
            covariance_col: Union[str,Integral] = "covariance",
            **kwargs: dict
        ) -> 'DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix':
        """
        Read a ``DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix`` from a ``pandas.DataFrame``.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.

        female2_col : str, Integral, default = "female2"
            Name or index of the column from which to read female 2 taxa names.

        female2_grp_col : str, Integral, None, default = "female2_grp"
            Name or index of the column from which to read female 2 taxa groups.

        male2_col : str, Integral, None, default = "male2"
            Name or index of the column from which to read male 2 taxa names.

        male2_grp_col : str, Integral, None, default = "male2_grp"
            Name or index of the column from which to read male 2 taxa names.

        female1_col : str, Integral, default = "female1"
            Name or index of the column from which to read female 1 taxa names.

        female1_grp_col : str, Integral, None, default = "female1_grp"
            Name or index of the column from which to read female 1 taxa groups.

        male1_col : str, Integral, None, default = "male1"
            Name or index of the column from which to read male 1 taxa names.

        male1_grp_col : str, Integral, None, default = "male1_grp"
            Name or index of the column from which to read male 1 taxa names.

        trait1_col : str, Integral, default = "trait"
            Name or index of the column from which to read trait taxa names.

        covariance_col : str, Integral, default = "covariance"
            Name or index of the column from which to read covariance taxa names.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            ``pandas.DataFrame``.
        
        Returns
        -------
        out : DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix
            A ``DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix`` read from a ``pandas.DataFrame``.
        """
        ### type checks and get df indices

        # df
        check_is_pandas_DataFrame(df, "df")

        # female2_col
        if isinstance(female2_col, str):
            check_pandas_DataFrame_has_column(df, "df", female2_col)
            female2_colix = df.columns.get_loc(female2_col)
        elif isinstance(female2_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", female2_col)
            female2_colix = female2_col
        else:
            check_is_str_or_Integral(female2_col, "female2_col")
        
        # female2_grp_col
        if female2_grp_col is not None:
            if isinstance(female2_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", female2_grp_col)
                female2_grp_colix = df.columns.get_loc(female2_grp_col)
            elif isinstance(female2_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", female2_grp_col)
                female2_grp_colix = female2_grp_col
            else:
                check_is_str_or_Integral(female2_grp_col, "female2_grp_col")
        
        # male2_col
        if isinstance(male2_col, str):
            check_pandas_DataFrame_has_column(df, "df", male2_col)
            male2_colix = df.columns.get_loc(male2_col)
        elif isinstance(male2_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", male2_col)
            male2_colix = male2_col
        else:
            check_is_str_or_Integral(male2_col, "male2_col")
        
        # male2_grp_col
        if male2_grp_col is not None:
            if isinstance(male2_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", male2_grp_col)
                male2_grp_colix = df.columns.get_loc(male2_grp_col)
            elif isinstance(male2_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", male2_grp_col)
                male2_grp_colix = male2_grp_col
            else:
                check_is_str_or_Integral(male2_grp_col, "male2_grp_col")
        
        # female1_col
        if isinstance(female1_col, str):
            check_pandas_DataFrame_has_column(df, "df", female1_col)
            female1_colix = df.columns.get_loc(female1_col)
        elif isinstance(female1_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", female1_col)
            female1_colix = female1_col
        else:
            check_is_str_or_Integral(female1_col, "female1_col")
        
        # female1_grp_col
        if female1_grp_col is not None:
            if isinstance(female1_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", female1_grp_col)
                female1_grp_colix = df.columns.get_loc(female1_grp_col)
            elif isinstance(female1_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", female1_grp_col)
                female1_grp_colix = female1_grp_col
            else:
                check_is_str_or_Integral(female1_grp_col, "female1_grp_col")
        
        # male1_col
        if isinstance(male1_col, str):
            check_pandas_DataFrame_has_column(df, "df", male1_col)
            male1_colix = df.columns.get_loc(male1_col)
        elif isinstance(male1_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", male1_col)
            male1_colix = male1_col
        else:
            check_is_str_or_Integral(male1_col, "male1_col")
        
        # male1_grp_col
        if male1_grp_col is not None:
            if isinstance(male1_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", male1_grp_col)
                male1_grp_colix = df.columns.get_loc(male1_grp_col)
            elif isinstance(male1_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", male1_grp_col)
                male1_grp_colix = male1_grp_col
            else:
                check_is_str_or_Integral(male1_grp_col, "male1_grp_col")
        
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
        female2_data  = df.iloc[:,female2_colix ].to_numpy(dtype = object)
        male2_data    = df.iloc[:,male2_colix   ].to_numpy(dtype = object)
        female1_data  = df.iloc[:,female1_colix ].to_numpy(dtype = object)
        male1_data    = df.iloc[:,male1_colix   ].to_numpy(dtype = object)
        trait1_data    = df.iloc[:,trait1_colix   ].to_numpy(dtype = object)
        trait2_data    = df.iloc[:,trait2_colix   ].to_numpy(dtype = object)
        covariance_data = df.iloc[:,covariance_colix].to_numpy(dtype = float)

        # get unique female2, male2, female1, male1 taxa (type numpy.ndarray)
        female2_taxa, female2_taxaix = numpy.unique(female2_data, return_index = True)
        male2_taxa, male2_taxaix     = numpy.unique(male2_data, return_index = True)
        female1_taxa, female1_taxaix = numpy.unique(female1_data, return_index = True)
        male1_taxa, male1_taxaix     = numpy.unique(male1_data, return_index = True)

        # combine female2, male2, female1, male1 taxa names (type numpy.ndarray)
        taxa = numpy.union1d(
            numpy.union1d(female2_taxa, male2_taxa),
            numpy.union1d(female1_taxa, male1_taxa)
        )

        # allocate female2, male2, female1, male1 indices
        female2ix = numpy.empty(len(female2_data), dtype = int)
        male2ix   = numpy.empty(len(male2_data  ), dtype = int)
        female1ix = numpy.empty(len(female1_data), dtype = int)
        male1ix   = numpy.empty(len(male1_data  ), dtype = int)

        # calculate female2, male2, female1, male1 indices
        for i,taxon in enumerate(taxa):
            female2ix[female2_data == taxon] = i
            male2ix  [  male2_data == taxon] = i
            female1ix[female1_data == taxon] = i
            male1ix  [  male1_data == taxon] = i
        
        # calculate unique trait values, trait indices
        trait1, trait1ix = numpy.unique(trait1_data, return_inverse = True)
        trait2, trait2ix = numpy.unique(trait2_data, return_inverse = True)

        trait = numpy.union1d(trait1, trait2)

        # get optional taxa group data
        taxa_grp = None
        if female2_grp_col is not None:
            if taxa_grp is None:
                taxa_grp = numpy.empty(len(taxa), dtype = int)
            female2_grp_data = df.iloc[:,female2_grp_colix].to_numpy(dtype = int)
            for i,taxon in enumerate(taxa):
                if taxon in female2_taxa:
                    ix = female2_taxaix[female2_taxa == taxon][0]
                    taxa_grp[i] = female2_grp_data[ix]
        if male2_grp_col is not None:
            if taxa_grp is None:
                taxa_grp = numpy.empty(len(taxa), dtype = int)
            male2_grp_data = df.iloc[:,male2_grp_colix].to_numpy(dtype = int)
            for i,taxon in enumerate(taxa):
                if taxon in male2_taxa:
                    ix = male2_taxaix[male2_taxa == taxon][0]
                    taxa_grp[i] = male2_grp_data[ix]
        if female1_grp_col is not None:
            if taxa_grp is None:
                taxa_grp = numpy.empty(len(taxa), dtype = int)
            female1_grp_data = df.iloc[:,female1_grp_colix].to_numpy(dtype = int)
            for i,taxon in enumerate(taxa):
                if taxon in female1_taxa:
                    ix = female1_taxaix[female1_taxa == taxon][0]
                    taxa_grp[i] = female1_grp_data[ix]
        if male1_grp_col is not None:
            if taxa_grp is None:
                taxa_grp = numpy.empty(len(taxa), dtype = int)
            male1_grp_data = df.iloc[:,male1_grp_colix].to_numpy(dtype = int)
            for i,taxon in enumerate(taxa):
                if taxon in male1_taxa:
                    ix = male1_taxaix[male1_taxa == taxon][0]
                    taxa_grp[i] = male1_grp_data[ix]

        # get array dimensions
        nfemale2 = len(taxa)
        nmale2 = len(taxa)
        nfemale1 = len(taxa)
        nmale1 = len(taxa)
        ntrait = len(trait)

        # allocate NaN array for covariance matrix
        mat = numpy.full((nfemale2,nmale2,nfemale1,nmale1,ntrait), numpy.nan, dtype = float)

        # overwrite NaN values with covariance values
        mat[female2ix,male2ix,female1ix,male1ix,trait1ix,trait2ix] = covariance_data

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
            female2_col: Union[str,Integral] = "female2",
            female2_grp_col: Optional[Union[str,Integral]] = "female2_grp",
            male2_col: Union[str,Integral] = "male2",
            male2_grp_col: Optional[Union[str,Integral]] = "male2_grp",
            female1_col: Union[str,Integral] = "female1",
            female1_grp_col: Optional[Union[str,Integral]] = "female1_grp",
            male1_col: Union[str,Integral] = "male1",
            male1_grp_col: Optional[Union[str,Integral]] = "male1_grp",
            trait1_col: Union[str,Integral] = "trait1",
            trait2_col: Union[str,Integral] = "trait2",
            covariance_col: Union[str,Integral] = "covariance",
            sep: str = ',',
            header: int = 0,
            **kwargs: dict
        ) -> 'DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix':
        """
        Read an object from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.

        female2_col : str, Integral, default = "female2"
            Name or index of the column from which to read female 2 taxa names.

        female2_grp_col : str, Integral, None, default = "female2_grp"
            Name or index of the column from which to read female 2 taxa groups.

        male2_col : str, Integral, default = "male2"
            Name or index of the column from which to read male 2 taxa names.

        male2_grp_col : str, Integral, None, default = "male2_grp"
            Name or index of the column from which to read male 2 taxa groups.

        female1_col : str, Integral, default = "female1"
            Name or index of the column from which to read female 1 taxa names.

        female1_grp_col : str, Integral, None, default = "female1_grp"
            Name or index of the column from which to read female 1 taxa groups.

        male1_col : str, Integral, None, default = "male1"
            Name or index of the column from which to read male 1 taxa names.

        male1_grp_col : str, Integral, None, default = "male1"
            Name or index of the column from which to read male 1 taxa names.

        trait1_col : str, Integral, default = "trait"
            Name or index of the column from which to read trait taxa names.

        covariance_col : str, Integral, default = "covariance"
            Name or index of the column from which to read covariance taxa names.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix
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
            female2_col = female2_col,
            female2_grp_col = female2_grp_col,
            male2_col = male2_col,
            male2_grp_col = male2_grp_col,
            female1_col = female1_col,
            female1_grp_col = female1_grp_col,
            male1_col = male1_col,
            male1_grp_col = male1_grp_col,
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
        ) -> 'DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix':
        """
        Read ``DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix`` from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which ``DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix`` data is stored.
            If ``None``, ``DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix`` is read from base HDF5 group.

        Returns
        -------
        out : DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix
            A ``DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix`` read from file.
        """
        # call super function
        return super(DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix, cls).from_hdf5(
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
        ) -> 'DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix':
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
        ) -> 'DenseFourWayDHAdditiveProgenyGenicCovarianceMatrix':
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
        epgc = (0.25,0.25,0.25,0.25)            # parental contributions

        # gather allele frequencies within each taxon
        tafreq = pgmat.tafreq()                 # (n,p) allele frequencies within taxon
        u = algmod.u_a                          # (p,t) marker effect coefficients
        
        # allocate a square matrix for each pairwise covariance
        cov_a = numpy.empty(
            (ntaxa,ntaxa,ntrait,ntrait),               # (n,n,t) covariance matrix
            dtype = float
        )

        # for each mate pair (including selfs)
        for female2 in range(0,ntaxa):
            for male2 in range(0,ntaxa):
                for female1 in range(0,ntaxa):
                    for male1 in range(0,female1):
                        # calculate the cross allele frequency
                        # (n,p)[(4,),:] -> (4,p)
                        # (4,) . (4,p) -> (p,)
                        p = numpy.dot(epgc, tafreq[(female2,male2,female1,male1),:])

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
                        cov_a[female2,male2,female1,male1,:,:] = cov
                        cov_a[female2,male2,male1,female1,:,:] = cov

        # construct output
        out = cls(
            mat = cov_a,
            taxa = pgmat.taxa,
            taxa_grp = pgmat.taxa_grp,
            trait = algmod.trait
        )

        return out
