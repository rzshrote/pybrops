import numpy

from pybropt.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol
from pybropt.core.error import check_is_str
from pybropt.core.error import check_is_array_like
from pybropt.popgen.ptdf import check_is_PhenotypeDataFrame
from pybropt.popgen.gmat import check_is_GenotypeMatrix
from pybropt.popgen.bvmat import DenseEstimatedBreedingValueMatrix

class MeanPhenotypicBreedingValue(BreedingValueProtocol):
    """docstring for MeanPhenotypicBreedingValue."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, taxa_col, trait_col, **kwargs):
        """
        Constructor for the concrete class Mean MeanPhenotypicBreedingValue.

        Calculate breeding values by taking the mean performance.

        Parameters
        ----------
        taxa_col : str
            Name of taxa column.
        trait_col : str, array_like
            Trait column names
        """
        super(MeanPhenotypicBreedingValue, self).__init__(**kwargs)

        # check taxa_col
        check_is_str(taxa_col, "taxa_col")
        self.taxa_col = taxa_col

        # check trait_col
        if isinstance(trait_col, str):  # if string
            trait_col = [trait_col]     # convert to list with single string
        check_is_array_like(trait_col, "trait_col")
        self.trait_col = trait_col

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def estimate(self, ptobj, gtobj, miscout = None, **kwargs):
        """
        Estimate breeding values by taking the mean performance.

        Parameters
        ----------
        ptobj : PhenotypeDataFrame
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
        check_is_PhenotypeDataFrame(ptobj, "ptobj")
        check_is_GenotypeMatrix(gtobj, "gtobj")

        # get taxa vector
        taxa_gt = gtobj.taxa
        taxa_pt = ptobj.col_data(name = self.taxa_col)

        # get unique phenotype taxa
        taxa_pt_uniq, taxa_pt_inv, taxa_pt_cnt = numpy.unique(
            taxa_pt,
            return_inverse = True,
            return_counts = True
        )

        # check to make sure we have all genotypes
        if not numpy.all(numpy.in1d(taxa_gt, taxa_pt_uniq)):
            raise ValueError("unable to align taxa: not all taxa present in gtobj present in ptobj")

        # get taxa selection index
        tselix = [numpy.flatnonzero(taxa_pt_uniq == e) for e in taxa_gt]
        tselix = numpy.concatenate(tselix)

        # taxa values
        mat = []

        # calculate means for each genotype
        for t in self.trait_col:
            tval = ptobj.col_data(name = t)[0]                  # get trait values for each taxa
            tsum = numpy.bincount(taxa_pt_inv, weights = tval)  # get trait sums for each taxa
            tmean = tsum / taxa_pt_cnt                          # get trait means for each taxa
            tmean = tmean[tselix]                               # select correct taxa
            mat.append(tmean)                                   # append trait means for selected taxa

        # stack results into an array
        mat = numpy.stack(mat, axis = 1)    # add each vector as a column

        # construct breeding value matrix
        bvmat = DenseEstimatedBreedingValueMatrix.from_numpy(
            a = mat,
            taxa = gtobj.taxa,
            taxa_grp = gtobj.taxa_grp,
            trait = numpy.object_(self.trait_col),
            **kwargs
        )

        return bvmat
