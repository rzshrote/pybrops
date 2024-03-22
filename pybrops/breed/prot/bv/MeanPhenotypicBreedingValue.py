"""
Module for estimating breeding values using the mean across all environments.
"""

from typing import Iterable
from typing import Optional
from typing import Union
import numpy
import pandas

from pybrops.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_column
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_columns
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.bvmat.DenseEstimatedBreedingValueMatrix import DenseEstimatedBreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.GenotypeMatrix import check_GenotypeMatrix_has_taxa
from pybrops.popgen.gmat.GenotypeMatrix import check_is_GenotypeMatrix

class MeanPhenotypicBreedingValue(BreedingValueProtocol):
    """
    Class implementing estimation of breeding values by taking the mean across
    all environments.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            taxa_col: str, 
            taxa_grp_col: Union[str,None],
            trait_cols: Union[str,Iterable], 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class Mean MeanPhenotypicBreedingValue.

        Calculate breeding values by taking the mean performance.

        Parameters
        ----------
        taxa_col : str
            Name of taxa column and/or names along which to group phenotypes along.
        trait_col : str, Iterable
            Name of trait column(s) along which to take the mean.
        kwargs : dict
            Additional keyword arguments.
        """
        # assignments
        self.taxa_col = taxa_col
        self.taxa_grp_col = taxa_grp_col
        self.trait_cols = trait_cols

    ############################ Object Properties #############################
    @property
    def taxa_col(self) -> str:
        """taxa_col."""
        return self._taxa_col
    @taxa_col.setter
    def taxa_col(self, value: str) -> None:
        """Set taxa_col."""
        check_is_str(value, "taxa_col")
        self._taxa_col = value

    @property
    def taxa_grp_col(self) -> Union[str,None]:
        """taxa_grp_col."""
        return self._taxa_grp_col
    @taxa_grp_col.setter
    def taxa_grp_col(self, value: Union[str,None]) -> None:
        """Set taxa_grp_col."""
        if value is not None:
            check_is_str(value, "taxa_grp_col")
        self._taxa_grp_col = value

    @property
    def trait_cols(self) -> list:
        """trait_col."""
        return self._trait_cols
    @trait_cols.setter
    def trait_cols(self, value: Union[str,Iterable]) -> None:
        """Set trait_col."""
        if isinstance(value, str):
            value = [value]
        elif isinstance(value, Iterable):
            value = list(value)
        else:
            raise TypeError("variable 'trait_col' must be of type str or an Iterable")
        self._trait_cols = value

    ############################## Object Methods ##############################
    def estimate(
            self, 
            ptobj: pandas.DataFrame, 
            gtobj: Optional[GenotypeMatrix] = None, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> BreedingValueMatrix:
        """
        Estimate breeding values by taking the mean performance.

        Parameters
        ----------
        ptobj : pandas.DataFrame
            An object containing phenotype data. Must be a phenotype data frame.
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a genotype matrix.
            This genotype object is used to align taxa orders in the returned
            BreedingValueMatrix.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : BreedingValueMatrix
            A matrix of breeding values.
        """
        # check arguments
        check_is_pandas_DataFrame(ptobj, "ptobj")
        check_pandas_DataFrame_has_column(ptobj, "ptobj", self.taxa_col)
        check_pandas_DataFrame_has_columns(ptobj, "ptobj", *self.trait_cols)
        if self.taxa_grp_col is not None:
            check_pandas_DataFrame_has_column(ptobj, "ptobj", self.taxa_grp_col)

        # get the columns to groupby
        by = [self.taxa_col]
        if self.taxa_grp_col is not None:
            by.append(self.taxa_grp_col)

        # take the means within each group for each trait
        agg_df = ptobj.groupby(
            by, 
            as_index = False, # do not use ``by`` as index
        ).agg(
            # generate ``trait: "mean"`` pairs
            dict((trait,"mean") for trait in self.trait_cols)
        )

        # if a Genotype matrix is not provided, then create matrix as-is
        if gtobj is None:
            mat = agg_df[self.trait_cols].to_numpy(dtype=float)
            taxa = agg_df[self.taxa_col].to_numpy(dtype=object)
            taxa_grp = None if self.taxa_grp_col is None else agg_df[self.taxa_grp_col].to_numpy(dtype=int)
            trait = numpy.array(self.trait_cols, dtype = object)
            out = DenseEstimatedBreedingValueMatrix.from_numpy(
                mat = mat,
                taxa = taxa,
                taxa_grp = taxa_grp,
                trait = trait,
            )
            return out

        # type and value checks
        check_is_GenotypeMatrix(gtobj, "gtobj")
        check_GenotypeMatrix_has_taxa(gtobj, "gtobj")

        # get calculated values
        agg_df_mat = agg_df[self.trait_cols].to_numpy()
        agg_df_taxa = agg_df[self.taxa_col].to_numpy()

        # make a hash table to look up row indices
        agg_df_taxa_hashtable = dict(zip(agg_df_taxa,range(len(agg_df_taxa))))

        # get the number of taxa and traits
        ntaxa = gtobj.ntaxa
        ntrait = len(self.trait_cols)

        # allocate matrix with NaN values
        mat = numpy.full((ntaxa,ntrait), numpy.nan, dtype = float)

        # for each value in the GenotypeMatrix, get the indices
        for i,taxon in enumerate(gtobj.taxa):
            # lookup the index
            try:
                ix = agg_df_taxa_hashtable[taxon]   # get index from hash table
                mat[i,:] = agg_df_mat[ix,:]         # copy data to output matrix
            except KeyError:
                continue
        
        # construct output
        out = DenseEstimatedBreedingValueMatrix.from_numpy(
            mat = mat,
            taxa = gtobj.taxa,
            taxa_grp = gtobj.taxa_grp,
            trait = numpy.array(self.trait_cols, dtype = object)
        )

        # copy metadata
        out.taxa_grp_name = gtobj.taxa_grp_name
        out.taxa_grp_stix = gtobj.taxa_grp_stix
        out.taxa_grp_spix = gtobj.taxa_grp_spix
        out.taxa_grp_len  = gtobj.taxa_grp_len

        return out
