"""
Module for estimating breeding values using the mean across all environments.
"""

from typing import Iterable, Optional, Sequence, Union
import numpy
import pandas

from pybrops.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_type_python import check_is_array_like
from pybrops.core.error.error_type_python import check_is_str
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.bvmat.DenseEstimatedBreedingValueMatrix import DenseEstimatedBreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix, check_GenotypeMatrix_has_taxa, check_is_GenotypeMatrix

class MeanPhenotypicBreedingValue(BreedingValueProtocol):
    """
    Class implementing estimation of breeding values by taking the mean across
    all environments.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            taxa_col: str, 
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

        # take the means within each group for each trait
        mean_ls = ["mean"] * len(self.trait_cols)
        agg_df = ptobj.groupby(self.taxa_col, as_index=False).agg(dict(zip(self.trait_cols,mean_ls)))

        # if a Genotype matrix is not provided, then create matrix as-is
        if gtobj is None:
            out = DenseEstimatedBreedingValueMatrix.from_numpy(
                mat = agg_df[self.trait_cols].to_numpy(),
                taxa = agg_df[self.taxa_col].to_numpy(),
                taxa_grp = None,
                trait = numpy.array(self.trait_cols, dtype = object)
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

        return out
