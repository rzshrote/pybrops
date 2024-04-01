"""
Module implementing classes and associated error checking routines for matrices
storing dense additive genetic covariance estimates calculated using four-way DH
formulae.
"""

import math
from numbers import Integral
from numbers import Real
from pathlib import Path
from typing import Optional
from typing import Union
import numpy
import pandas
import h5py
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_type_python import check_is_Integral_or_None
from pybrops.core.error.error_type_python import check_is_Integral_or_inf
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_type_python import check_is_str_or_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_column
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_column_index
from pybrops.core.util.array import flattenix
from pybrops.core.util.subroutines import srange
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.AdditiveLinearGenomicModel import check_is_AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.gmod.GenomicModel import check_is_GenomicModel
from pybrops.model.vmat.util import cov_D1s
from pybrops.model.vmat.util import cov_D2s
from pybrops.model.pcvmat.DenseAdditiveProgenyGeneticCovarianceMatrix import DenseAdditiveProgenyGeneticCovarianceMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmap.GeneticMapFunction import check_is_GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import check_is_PhasedGenotypeMatrix

class DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix(
        DenseAdditiveProgenyGeneticCovarianceMatrix,
    ):
    """
    A concrete class for dense additive genetic covariance matrices calculated
    for four-way DH progenies.

    The purpose of this concrete class is to implement functionality for:
        1) Genetic covariance estimation for four-way DH progenies.
        2) I/O for four-way DH progeny covariance matrices.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix.

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
        super(DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix, self).__init__(
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
    @DenseAdditiveProgenyGeneticCovarianceMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set pointer to raw numpy.ndarray object."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim(value, "mat", 6) # (ntaxa,ntaxa,ntaxa,ntaxa,ntrait,ntrait)
        self._mat = value

    ############## Square Metadata Properties ##############
    @DenseAdditiveProgenyGeneticCovarianceMatrix.square_taxa_axes.getter
    def square_taxa_axes(self) -> tuple:
        """Axis indices for taxa axes that are square."""
        return (0,1,2,3) # (female2, male2, female1, male1)

    @DenseAdditiveProgenyGeneticCovarianceMatrix.square_trait_axes.getter
    def square_trait_axes(self) -> tuple:
        return (4,5) # (trait1, trait2) covariance matrix

    #################### Trait metadata ####################

    ######## Expected parental genome contributions ########
    @DenseAdditiveProgenyGeneticCovarianceMatrix.epgc.getter
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
        Export a DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix to a pandas.DataFrame.

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

        trait_col : str, default = "trait"
            Name of the column to which to write trait taxa names.

        variance_col : str, default = "covariance"
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
        check_is_str(covariance_col, "variance_col")

        # if taxa is None, give default names of TaxonXXX, where X is the number.
        if self.taxa is None:
            nzero = math.ceil(math.log10(self.ntaxa))+1
            taxa = numpy.array(["Taxon"+str(e).zfill(nzero) for e in range(self.ntaxa)], dtype = object)
        else:
            taxa = self.taxa

        # if trait is None, give default names of TraitXXX, where X is the number.
        if self.trait is None:
            nzero = math.ceil(math.log10(self.ntrait))+1
            trait = numpy.array(["Trait"+str(e).zfill(nzero) for e in range(self.ntrait)], dtype = object)
        else:
            trait = self.trait

        # calculate flattened array and corresponding axis indices
        flatmat, (female2ix,male2ix,female1ix,male1ix,trait1ix,trait2ix) = flattenix(self.mat)

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
        Write a DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix to a CSV file.

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

        trait_col : str, default = "trait"
            Name of the column to which to write trait taxa names.

        variance_col : str, default = "covariance"
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
        # convert DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix to pandas.DataFrame
        df = self.to_pandas(
            female2_col     = female2_col,
            female2_grp_col = female2_grp_col,
            male2_col       = male2_col,
            male2_grp_col   = male2_grp_col,
            female1_col     = female1_col,
            female1_grp_col = female1_grp_col,
            male1_col       = male1_col,
            male1_grp_col   = male1_grp_col,
            trait1_col      = trait1_col,
            trait2_col      = trait2_col,
            covariance_col  = covariance_col,
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
        Write ``DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix`` to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            HDF5 file name or HDF5 file stream to which to write.

        groupname : str, None
            If ``str``, an HDF5 group name under which the ``DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix`` data is stored.
            If ``None``, the ``DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix`` is written to the base HDF5 group.
        
        overwrite : bool
            Whether to overwrite data fields if they are present in the HDF5 file.
        """
        # call super function
        super(DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix, self).to_hdf5(
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
        ) -> 'DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix':
        """
        Read a ``DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix`` from a ``pandas.DataFrame``.

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

        trait_col : str, Integral, default = "trait"
            Name or index of the column from which to read trait taxa names.

        variance_col : str, Integral, default = "covariance"
            Name or index of the column from which to read covariance taxa names.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            ``pandas.DataFrame``.
        
        Returns
        -------
        out : DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix
            A ``DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix`` read from a ``pandas.DataFrame``.
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
        
        # variance_col
        if isinstance(covariance_col, str):
            check_pandas_DataFrame_has_column(df, "df", covariance_col)
            covariance_colix = df.columns.get_loc(covariance_col)
        elif isinstance(covariance_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", covariance_col)
            covariance_colix = covariance_col
        else:
            check_is_str_or_Integral(covariance_col, "variance_col")
        
        # get required data columns (type numpy.ndarray)
        female2_data    = df.iloc[:,female2_colix   ].to_numpy(dtype = object)
        male2_data      = df.iloc[:,male2_colix     ].to_numpy(dtype = object)
        female1_data    = df.iloc[:,female1_colix   ].to_numpy(dtype = object)
        male1_data      = df.iloc[:,male1_colix     ].to_numpy(dtype = object)
        trait1_data     = df.iloc[:,trait1_colix    ].to_numpy(dtype = object)
        trait2_data     = df.iloc[:,trait2_colix    ].to_numpy(dtype = object)
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

        # get all unique traits
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
        mat = numpy.full((nfemale2,nmale2,nfemale1,nmale1,ntrait,ntrait), numpy.nan, dtype = float)

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
        ) -> 'DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix':
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

        trait_col : str, Integral, default = "trait"
            Name or index of the column from which to read trait taxa names.

        variance_col : str, Integral, default = "covariance"
            Name or index of the column from which to read covariance taxa names.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix
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
        ) -> 'DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix':
        """
        Read ``DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix`` from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which ``DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix`` data is stored.
            If ``None``, ``DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix`` is read from base HDF5 group.

        Returns
        -------
        out : DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix
            A ``DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix`` read from file.
        """
        # call super function
        return super(DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix, cls).from_hdf5(
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
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Union[Integral,Real],
            gmapfn: GeneticMapFunction,
            **kwargs: dict
        ) -> 'DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix':
        """
        Calculate a symmetrical tensor of progeny variances for each possible
        3-way cross between *inbred* individuals.

        Parameters
        ----------
        gmod : GenomicModel
            Genomic Model with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : Integral
            Number of cross patterns to simulate for genetic covariance
            estimation.
        nprogeny : Integral
            Number of progeny to simulate per cross to estimate genetic
            covariance.
        nself : Integral, Real
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.

            +-----------------+-------------------------+
            | Example         | Description             |
            +=================+=========================+
            | ``nself = 0``   | Derive gametes from F1  |
            +-----------------+-------------------------+
            | ``nself = 1``   | Derive gametes from F2  |
            +-----------------+-------------------------+
            | ``nself = 2``   | Derive gametes from F3  |
            +-----------------+-------------------------+
            | ``...``         | etc.                    |
            +-----------------+-------------------------+
            | ``nself = inf`` | Derive gametes from SSD |
            +-----------------+-------------------------+
        gmapfn : GeneticMapFunction
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix
            A matrix of additive genetic covariance estimations.
        """
        # type checks
        check_is_GenomicModel(gmod, "gmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral_or_inf(nself, "nself")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")

        # if genomic model is an additive linear genomic model, then use specialized routine
        if isinstance(gmod, AdditiveLinearGenomicModel):
            return cls.from_algmod(
                algmod = gmod, 
                pgmat = pgmat, 
                ncross = ncross, 
                nprogeny = nprogeny, 
                nself = nself, 
                gmapfn = gmapfn, 
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
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Union[Integral,Real], 
            gmapfn: GeneticMapFunction, 
            mem: Union[Integral,None] = 1024
        ) -> 'DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix':
        """
        Calculate a symmetrical tensor of progeny variances for each possible
        4-way cross between *inbred* individuals.
        Calculations are derived from Lehermeier et al. (2019).

        Parameters
        ----------
        algmod : AdditiveLinearGenomicModel
            AdditiveLinearGenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : Integral
            Number of cross patterns to simulate for genetic covariance
            estimation.
        nprogeny : Integral
            Number of progeny to simulate per cross to estimate genetic
            covariance.
        nself : Integral
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.

            +-----------------+-------------------------+
            | Example         | Description             |
            +=================+=========================+
            | ``nself = 0``   | Derive gametes from F1  |
            +-----------------+-------------------------+
            | ``nself = 1``   | Derive gametes from F2  |
            +-----------------+-------------------------+
            | ``nself = 2``   | Derive gametes from F3  |
            +-----------------+-------------------------+
            | ``...``         | etc.                    |
            +-----------------+-------------------------+
            | ``nself = inf`` | Derive gametes from SSD |
            +-----------------+-------------------------+
        gmapfn : GeneticMapFunction
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        mem : Integral, default = 1024
            Memory chunk size to use during matrix operations. If ``None``,
            then memory chunk size is not limited.

            WARNING: Setting ``mem = None`` might result in memory allocation
            errors! For reference, ``mem = 1024`` refers to a matrix of size
            1024x1024, which needs about 8.5 MB of storage. Matrices of course
            need a quadratic amount of memory: :math:`O(n^2)`.

        Returns
        -------
        out : ProgenyGeneticCovarianceMatrix
            A matrix of additive genetic covariance estimations.
        """
        # type checks
        check_is_AdditiveLinearGenomicModel(algmod, "algmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_Integral(ncross, "ncross")
        check_is_Integral(nprogeny, "nprogeny")
        check_is_Integral_or_inf(nself, "nself")
        check_is_GeneticMapFunction(gmapfn, "gmapfn")
        check_is_Integral_or_None(mem, "mem")
        
        # check for chromosome grouping
        if not pgmat.is_grouped_vrnt():
            raise ValueError("pgmat must be grouped along the vrnt axis")

        # check for genetic positions
        if pgmat.vrnt_genpos is None:
            raise ValueError("pgmat must have genetic positions")

        # gather shapes of data input
        ntrait = algmod.ntrait                  # number ot traits (t)
        ntaxa = pgmat.ntaxa                     # number of individuals (n)

        # gather data pointers
        geno = pgmat.mat                        # (m,n,p) genotype matrix
        chrgrp_stix = pgmat.vrnt_chrgrp_stix    # (g,) chromosome group start indices
        chrgrp_spix = pgmat.vrnt_chrgrp_spix    # (g,) chromosome group stop indices
        genpos = pgmat.vrnt_genpos              # (p,) marker genetic positions
        u = algmod.u_a                          # (p,t) marker effect coefficients

        # allocate a square matrix for each pairwise covariance
        varA_mat = numpy.zeros(
            (ntaxa, ntaxa, ntaxa, ntaxa, ntrait),   # (n,n,n,n,t) covariance matrix
            dtype = float
        )

        # for each linkage group
        for lst, lsp in zip(chrgrp_stix, chrgrp_spix):
            # determine memory chunk step (for iterators)
            step = (lsp - lst) if mem is None else mem
            # for each computational chunk
            for rst,rsp in zip(range(lst,lsp,step),srange(lst+step,lsp,step)):
                for cst,csp in zip(range(lst,lsp,step),srange(lst+step,lsp,step)):
                    # create sparse meshgrid indicating where genetic positions are
                    gi, gj = numpy.meshgrid(
                        genpos[rst:rsp],        # row block genetic positions
                        genpos[cst:csp],        # column block genetic positions
                        indexing = 'ij',        # use ij indexing
                        sparse = True           # generate a spare array tuple for speed
                    )

                    # calculate recombination probability matrix for chunk
                    r = gmapfn.mapfn(numpy.abs(gi - gj)) # (rb,cb) covariance matrix

                    # calculate a D1 matrix; this is specific to mating scheme
                    D1 = cov_D1s(r, nself)  # (rb,cb) covariance matrix

                    # calculate a D2 matrix; this is specific to mating scheme
                    D2 = cov_D2s(r, nself)  # (rb,cb) covariance matrix

                    # get marker coefficients for rows and columns
                    ru = u[rst:rsp].T # (rb,t)' -> (t,rb)
                    cu = u[cst:csp].T # (cb,t)' -> (t,cb)

                    # for each 4-way cross (excluding selfs)
                    # subscript codes:
                    #   1 = female 2
                    #   2 = male 2
                    #   3 = female 1
                    #   4 = male 1
                    # TODO: make this more computationally efficient.
                    #       operations can be simplified if you mirror sections
                    #       of the 4d matrix between blocks.
                    #       how to do this is difficult to engineer
                    for female2 in range(0,ntaxa):              # varA block index (change me for efficiency?)
                        for male2 in range(0,ntaxa):            # varA slice index (change me for efficiency?)
                            # calculate genotype differences for row, col
                            rdgeno21 = geno[0,male2,rst:rsp] - geno[0,female2,rst:rsp] # (rb,)
                            cdgeno21 = geno[0,male2,cst:csp] - geno[0,female2,cst:csp] # (cb,)

                            # calculate effect differences
                            reffect21 = rdgeno21 * ru # (rb,)*(t,rb) -> (t,rb)
                            ceffect21 = cdgeno21 * cu # (cb,)*(t,cb) -> (t,cb)

                            # calculate the dot product for each trait to get a
                            # partial covariance sum for male2-female2
                            # (t,rb)@(rb,cb)@(cb,t) -> (t,t)
                            varA_part21 = reffect21 @ D2 @ ceffect21.T

                            for female1 in range(0,ntaxa):      # varA row index
                                # calculate genotype differences for row, col
                                rdgeno31 = geno[0,female1,rst:rsp] - geno[0,female2,rst:rsp] # (rb,)
                                cdgeno31 = geno[0,female1,cst:csp] - geno[0,female2,cst:csp] # (cb,)
                                rdgeno32 = geno[0,female1,rst:rsp] - geno[0,male2,rst:rsp] # (rb,)
                                cdgeno32 = geno[0,female1,cst:csp] - geno[0,male2,cst:csp] # (cb,)

                                # calculate effect differences
                                reffect31 = rdgeno31 * ru # (rb,)*(t,rb) -> (t,rb)
                                ceffect31 = cdgeno31 * cu # (cb,)*(t,cb) -> (t,cb)
                                reffect32 = rdgeno32 * ru # (rb,)*(t,rb) -> (t,rb)
                                ceffect32 = cdgeno32 * cu # (cb,)*(t,cb) -> (t,cb)

                                # calculate varA part for female1
                                # (t,rb)@(rb,cb)@(cb,t) -> (t,cb)
                                varA_part31 = reffect31 @ D1 @ ceffect31.T
                                varA_part32 = reffect32 @ D1 @ ceffect32.T

                                # only do lower triangle since symmetrical within each slice
                                for male1 in range(0,female1):  # varA col index
                                    # calculate genotype differences for row, col
                                    rdgeno41 = geno[0,male1,rst:rsp] - geno[0,female2,rst:rsp] # (rb,)
                                    cdgeno41 = geno[0,male1,cst:csp] - geno[0,female2,cst:csp] # (cb,)
                                    rdgeno42 = geno[0,male1,rst:rsp] - geno[0,male2,rst:rsp] # (rb,)
                                    cdgeno42 = geno[0,male1,cst:csp] - geno[0,male2,cst:csp] # (cb,)
                                    rdgeno43 = geno[0,male1,rst:rsp] - geno[0,female1,rst:rsp] # (rb,)
                                    cdgeno43 = geno[0,male1,cst:csp] - geno[0,female1,cst:csp] # (cb,)

                                    # calculate effect differences
                                    reffect41 = rdgeno41 * ru # (rb,)*(t,rb) -> (t,rb)
                                    ceffect41 = cdgeno41 * cu # (cb,)*(t,cb) -> (t,cb)
                                    reffect42 = rdgeno42 * ru # (rb,)*(t,rb) -> (t,rb)
                                    ceffect42 = cdgeno42 * cu # (cb,)*(t,cb) -> (t,cb)
                                    reffect43 = rdgeno43 * ru # (rb,)*(t,rb) -> (t,rb)
                                    ceffect43 = cdgeno43 * cu # (cb,)*(t,cb) -> (t,cb)

                                    # calculate varA parts for crosses with male
                                    # (t,rb)@(rb,cb)@(cb,t) -> (t,cb)
                                    varA_part41 = reffect41 @ D1 @ ceffect41.T
                                    varA_part42 = reffect42 @ D1 @ ceffect42.T
                                    varA_part43 = reffect43 @ D2 @ ceffect43.T

                                    # calculate varA part for this matrix chunk
                                    # (t,t) + (t,t) + (t,t) + (t,t) + (t,t) + (t,t) -> (t,t)
                                    varA_part = varA_part21 + varA_part31 + varA_part32 + varA_part41 + varA_part42 + varA_part43

                                    # add this partial covariance to the lower triangle
                                    # (1,1,1,1,t,t) + (t,t) -> (1,1,1,1,t,t)
                                    varA_mat[female2,male2,female1,male1,:,:] += varA_part

        # divide entire matrix by 4 to get covariance per the equation
        # multiplication is faster computationally
        varA_mat *= 0.25

        # each matrix is symmetrical within a slice because exchanging female
        # and male orders is mathematically equivalent.
        # copy lower triangle to the upper since varA matrix is symmetrical within each slice
        for female1 in range(1, ntaxa):
            for male1 in range(0, female1):
                varA_mat[:,:,male1,female1,:,:] = varA_mat[:,:,female1,male1,:,:]

        # construct output
        out = cls(
            mat = varA_mat,
            taxa = pgmat.taxa,
            taxa_grp = pgmat.taxa_grp
        )

        return out



################################## Utilities ###################################
def check_is_DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix):
        raise TypeError("'{0}' must be a DenseFourWayDHAdditiveProgenyGeneticCovarianceMatrix".format(vname))
