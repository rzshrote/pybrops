import pandas

from . import PhenotypingProtocol
from pybropt.popgen.ptdf import PandasPhenotypeDataFrame

class TruePhenotyping(PhenotypingProtocol):
    """docstring for TruePhenotyping."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(TruePhenotyping, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def phenotype(self, pgmat, gpmod, **kwargs):
        """
        Phenotype a set of genotypes using a genomic prediction model.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes of the individuals to phenotype.
        gpmod : GenomicModel
            Genomic prediction model to use to determine phenotypes.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhenotypeDataFrame
            DataFrame containing phenotypes.
        """
        # gather true breeding values
        bvmat = gpmod.predict(pgmat)

        # gather pointers to raw matrices
        mat = bvmat.mat             # breeding values
        taxa = bvmat.taxa           # taxa names
        taxa_grp = bvmat.taxa_grp   # taxa groups
        trait = bvmat.trait         # trait names
        ndim = bvmat.mat_ndim       # number of matrix dimensions
        taxis = bvmat.trait_axis    # trait axis

        # perform error checks
        if taxa is None:
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' does not have taxa names")
        if any(e is None for e in taxa):
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' has taxa name(s) which are 'None'")
        if trait is None:
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' does not have trait names")
        if any(e is None for e in trait):
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' has trait name(s) which are 'None'")

        # construct dictionary
        data_dict = {}

        # add taxa names
        data_dict["taxa"] = taxa

        # if there are taxa groups, add group information
        if taxa_grp is not None:
            data_dict["taxa_grp"] = taxa_grp

        # add each trait and corresponding data
        for i,e in enumerate(trait):
            # construct matrix slice selection tuple
            t = tuple(i if a == taxis else slice(None) for a in range(ndim))
            # select data and put into dictionary
            data_dict[e] = mat[t]

        # construct pandas.DataFrame
        df = pandas.DataFrame(data_dict)

        # construct PandasPhenotypeDataFrame
        ptdf = PandasPhenotypeDataFrame(
            df = df,
            **kwargs
        )

        return ptdf
