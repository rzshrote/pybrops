from . import PhenotypingProtocol
from pybropt.popgen.ptdf import DictPhenotypeDataFrame

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
        bvmat = gpmod.gebv(pgmat)

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

        # construct data dictionary
        data_dict = {"taxa": taxa}

        # construct column analysis type dictionary
        col_analysis_type_dict = {"taxa": "factor(str)"}

        # construct column effect type dictionary
        col_analysis_effect_dict = {"taxa": "fixed"}

        # if there are taxa groups, add group information
        if taxa_grp is not None:
            data_dict.update(taxa_grp = taxa_grp)
            col_analysis_type_dict.update(taxa_grp = "factor(int)")
            col_analysis_effect_dict.update(taxa_grp = "fixed")

        # trait type
        if mat.dtype == "float64":
            tatype = "double"
        elif mat.dtype == "float32":
            mat = mat.astype("float64")
            tatype = "double"
        elif mat.dtype == "int8":
            mat = mat.astype("int32")
            tatype = "int"
        elif mat.dtype == "int16":
            mat = mat.astype("int32")
            tatype = "int"
        elif mat.dtype == "int32":
            tatype = "int"
        elif mat.dtype == "int64":
            mat = mat.astype("int32")
            tatype = "int"
        else:
            raise TypeError("unsupported breeding value data type")

        # add each trait and corresponding data
        for i,e in enumerate(trait):
            # construct matrix slice selection tuple
            t = tuple(i if a == taxis else slice(None) for a in range(ndim))
            data_dict[e] = mat[t].copy()                # select data and put into dictionary
            col_analysis_type_dict[e] = tatype          # add column analysis type
            col_analysis_effect_dict[e] = "response"    # add column analysis effect type

        # construct DictPhenotypeDataFrame
        ptdf = DictPhenotypeDataFrame(
            data = data_dict,
            col_grp = None,
            col_analysis_type = col_analysis_type_dict,
            col_analysis_effect = col_analysis_effect_dict,
            row_name = None,
            **kwargs
        )

        return ptdf
