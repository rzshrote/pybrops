"""
Module implementing a dense matrix with two square taxa axes and one trait axis
and associated error checking routines.
"""

import copy
import math
from numbers import Integral
from typing import Optional, Union
import numpy
import pandas
from pandas.core.api import DataFrame as DataFrame
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_type_python import check_is_str, check_is_str_or_Integral
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_column, check_pandas_DataFrame_has_column_index
from pybrops.core.io.PandasInputOutput import PandasInputOutput
from pybrops.core.mat.DenseSquareTaxaTraitMatrix import DenseSquareTaxaTraitMatrix
from pybrops.core.util.array import flattenix

class DenseSquare2TaxaTraitMatrix(
        DenseSquareTaxaTraitMatrix,
        PandasInputOutput,
    ):
    """
    A concrete class for dense matrices with two square taxa axes and one trait 
    axis.
    """

    ########################## Special Object Methods ##########################

    ##################### Constructor ######################
    def __init__(
            self,
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseSquare2TaxaTraitMatrix.
        
        Parameters
        ----------
        mat : numpy.ndarray
            Matrix used to construct the object.
        taxa : numpy.ndarray
            Taxa names.
        taxa_grp : numpy.ndarray
            Taxa groupings.
        trait : numpy.ndarray
            Trait labels.
        kwargs : dict
            Additional keyword arguments.
        """
        # call DenseSquareTaxaTraitMatrix constructor
        super(DenseSquare2TaxaTraitMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseSquare2TaxaTraitMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseSquare2TaxaTraitMatrix
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
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
            memo: dict,
        ) -> 'DenseSquare2TaxaTraitMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseSquare2TaxaTraitMatrix
        """
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
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

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseSquare2TaxaTraitMatrix':
        """
        Make a shallow copy of the Matrix.

        Returns
        -------
        out : DenseSquare2TaxaTraitMatrix
            A shallow copy of the original DenseSquare2TaxaTraitMatrix.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseSquare2TaxaTraitMatrix':
        """
        Make a deep copy of the Matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseSquare2TaxaTraitMatrix
            A deep copy of the original DenseSquare2TaxaTraitMatrix.
        """
        return copy.deepcopy(self, memo)

    ############################## Object Methods ##############################

    def to_pandas(
            self,
            taxa1_col: str = "taxa1",
            taxa1_grp_col: Optional[str] = "taxa1_grp",
            taxa2_col: str = "taxa2",
            taxa2_grp_col: Optional[str] = "taxa2_grp",
            trait_col: str = "trait", 
            value_col: str = "value",
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Export a DenseSquare2TaxaTraitMatrix to a pandas.DataFrame.

        Parameters
        ----------
        taxa1_col : str, default = "taxa1"
            Name of the column to which to write female taxa names.

        taxa1_grp_col : str, None, default = "taxa1_grp"
            Name of the column to which to write female taxa groups.

        taxa2_col : str, default = "taxa2"
            Name of the column to which to write male taxa names.

        taxa2_grp_col : str, None, default = "taxa2_grp"
            Name of the column to which to write male taxa groups.

        trait_col : str, default = "trait"
            Name of the column to which to write trait taxa names.

        value_col : str, default = "value"
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
        check_is_str(taxa1_col, "taxa1_col")
        if taxa1_grp_col is not None:
            check_is_str(taxa1_grp_col, "taxa1_grp_col")
        check_is_str(taxa2_col, "taxa2_col")
        if taxa2_grp_col is not None:
            check_is_str(taxa2_grp_col, "taxa2_grp_col")
        check_is_str(trait_col, "trait_col")

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
        flatmat, (femaleix,maleix,traitix) = flattenix(self.mat)

        # make dictionary to store output columns in specific column ordering
        out_dict = {}
        out_dict.update({taxa1_col: taxa[femaleix]})
        if taxa1_grp_col is not None:
            values = None if self.taxa_grp is None else self.taxa_grp[femaleix]
            out_dict.update({taxa1_grp_col: values})
        out_dict.update({taxa2_col: taxa[maleix]})
        if taxa2_grp_col is not None:
            values = None if self.taxa_grp is None else self.taxa_grp[maleix]
            out_dict.update({taxa2_grp_col: values})
        out_dict.update({trait_col: trait[traitix]})
        out_dict.update({value_col: flatmat})

        # create a pandas DataFrame from the data
        out = pandas.DataFrame(out_dict, **kwargs)

        return out

    ############################## Class Methods ###############################

    ###################### Matrix I/O ######################
    @classmethod
    def from_pandas(
            cls, 
            df: pandas.DataFrame, 
            taxa1_col: Union[str,Integral] = "female",
            taxa1_grp_col: Optional[Union[str,Integral]] = "female_grp",
            taxa2_col: Union[str,Integral] = "male",
            taxa2_grp_col: Optional[Union[str,Integral]] = "male_grp",
            trait_col: Union[str,Integral] = "trait1",
            value_col: Union[str,Integral] = "covariance",
            **kwargs: dict
        ) -> 'DenseSquare2TaxaTraitMatrix':
        """
        Read a ``DenseSquare2TaxaTraitMatrix`` from a ``pandas.DataFrame``.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.

        taxa1_col : str, Integral, default = "female"
            Name or index of the column from which to read female taxa names.

        taxa1_grp_col : str, Integral, None, default = "female"
            Name or index of the column from which to read female taxa groups.

        taxa2_col : str, Integral, None, default = "male"
            Name or index of the column from which to read male taxa names.

        taxa2_grp_col : str, Integral, None, default = "male"
            Name or index of the column from which to read male taxa names.

        trait_col : str, Integral, default = "trait"
            Name or index of the column from which to read trait taxa names.

        value_col : str, Integral, default = "covariance"
            Name or index of the column from which to read covariance taxa names.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            ``pandas.DataFrame``.
        
        Returns
        -------
        out : DenseSquare2TaxaTraitMatrix
            A ``DenseSquare2TaxaTraitMatrix`` read from a ``pandas.DataFrame``.
        """
        ### type checks and get df indices

        # df
        check_is_pandas_DataFrame(df, "df")

        # taxa1_col
        if isinstance(taxa1_col, str):
            check_pandas_DataFrame_has_column(df, "df", taxa1_col)
            taxa1_colix = df.columns.get_loc(taxa1_col)
        elif isinstance(taxa1_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", taxa1_col)
            taxa1_colix = taxa1_col
        else:
            check_is_str_or_Integral(taxa1_col, "taxa1_col")
        
        # taxa1_grp_col
        if taxa1_grp_col is not None:
            if isinstance(taxa1_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", taxa1_grp_col)
                taxa1_grp_colix = df.columns.get_loc(taxa1_grp_col)
            elif isinstance(taxa1_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", taxa1_grp_col)
                taxa1_grp_colix = taxa1_grp_col
            else:
                check_is_str_or_Integral(taxa1_grp_col, "taxa1_grp_col")
        
        # taxa2_col
        if isinstance(taxa2_col, str):
            check_pandas_DataFrame_has_column(df, "df", taxa2_col)
            taxa2_colix = df.columns.get_loc(taxa2_col)
        elif isinstance(taxa2_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", taxa2_col)
            taxa2_colix = taxa2_col
        else:
            check_is_str_or_Integral(taxa2_col, "taxa2_col")
        
        # taxa2_grp_col
        if taxa2_grp_col is not None:
            if isinstance(taxa2_grp_col, str):
                check_pandas_DataFrame_has_column(df, "df", taxa2_grp_col)
                taxa2_grp_colix = df.columns.get_loc(taxa2_grp_col)
            elif isinstance(taxa2_grp_col, Integral):
                check_pandas_DataFrame_has_column_index(df, "df", taxa2_grp_col)
                taxa2_grp_colix = taxa2_grp_col
            else:
                check_is_str_or_Integral(taxa2_grp_col, "taxa2_grp_col")
        
        # trait_col
        if isinstance(trait_col, str):
            check_pandas_DataFrame_has_column(df, "df", trait_col)
            trait_colix = df.columns.get_loc(trait_col)
        elif isinstance(trait_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", trait_col)
            trait_colix = trait_col
        else:
            check_is_str_or_Integral(trait_col, "trait_col")
        
        # value_col
        if isinstance(value_col, str):
            check_pandas_DataFrame_has_column(df, "df", value_col)
            value_colix = df.columns.get_loc(value_col)
        elif isinstance(value_col, Integral):
            check_pandas_DataFrame_has_column_index(df, "df", value_col)
            value_colix = value_col
        else:
            check_is_str_or_Integral(value_col, "value_col")
        
        # get required data columns (type numpy.ndarray)
        female_data     = df.iloc[:,taxa1_colix].to_numpy(dtype = object)
        male_data       = df.iloc[:,taxa2_colix].to_numpy(dtype = object)
        trait_data      = df.iloc[:,trait_colix].to_numpy(dtype = object)
        covariance_data = df.iloc[:,value_colix].to_numpy(dtype = float)

        # get unique female, male taxa (type numpy.ndarray)
        female_taxa, female_taxaix = numpy.unique(female_data, return_index = True)
        male_taxa, male_taxaix     = numpy.unique(male_data, return_index = True)

        # combine female and male taxa names (type numpy.ndarray)
        taxa = numpy.union1d(female_taxa, male_taxa)

        # allocate female, male indices
        femaleix = numpy.empty(len(female_data), dtype = int)
        maleix   = numpy.empty(len(male_data  ), dtype = int)

        # calculate female, male indices
        for i,taxon in enumerate(taxa):
            femaleix[female_data == taxon] = i
            maleix[male_data == taxon] = i
        
        # calculate unique trait values, trait indices
        trait, traitix = numpy.unique(trait_data, return_inverse = True)

        # get optional taxa group data
        taxa_grp = None
        if taxa1_grp_col is not None:
            if taxa_grp is None:
                taxa_grp = numpy.empty(len(taxa), dtype = int)
            female_grp_data = df.iloc[:,taxa1_grp_colix].to_numpy(dtype = int)
            for i,taxon in enumerate(taxa):
                if taxon in female_taxa:
                    ix = female_taxaix[female_taxa == taxon][0]
                    taxa_grp[i] = female_grp_data[ix]
        if taxa2_grp_col is not None:
            if taxa_grp is None:
                taxa_grp = numpy.empty(len(taxa), dtype = int)
            male_grp_data = df.iloc[:,taxa2_grp_colix].to_numpy(dtype = int)
            for i,taxon in enumerate(taxa):
                if taxon in male_taxa:
                    ix = male_taxaix[male_taxa == taxon][0]
                    taxa_grp[i] = male_grp_data[ix]

        # get array dimensions
        nfemale = len(taxa)
        nmale = len(taxa)
        ntrait = len(trait)

        # allocate NaN array for covariance matrix
        mat = numpy.full((nfemale,nmale,ntrait,ntrait), numpy.nan, dtype = float)

        # overwrite NaN values with covariance values
        mat[femaleix,maleix,traitix] = covariance_data

        # construct an object
        out = cls(
            mat = mat, 
            taxa = taxa, 
            taxa_grp = taxa_grp, 
            trait = trait, 
        )

        return out


    ############################## Static Methods ##############################
