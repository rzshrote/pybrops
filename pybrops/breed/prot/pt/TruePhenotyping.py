import numpy

from pybrops.breed.prot.pt.PhenotypingProtocol import PhenotypingProtocol
from pybrops.popgen.ptdf.DictPhenotypeDataFrame import DictPhenotypeDataFrame
from pybrops.model.gmod.GenomicModel import check_is_GenomicModel
from pybrops.core.error import error_readonly
from pybrops.popgen.gmat.PhasedGenotypeMatrix import check_is_PhasedGenotypeMatrix

class TruePhenotyping(PhenotypingProtocol):
    """docstring for TruePhenotyping."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, gpmod, **kwargs):
        """
        Constructor for the concrete class TruePhenotyping.

        Parameters
        ----------
        gpmod : GenomicModel
            Genomic prediction model to use to determine phenotypes.
        kwargs : dict
            Additional keyword arguments
        """
        super(TruePhenotyping, self).__init__(**kwargs)
        self.gpmod = gpmod

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############### Genomic Model Properties ###############
    def gpmod():
        doc = "Genomic prediction model."
        def fget(self):
            """Get genomic prediction model"""
            return self._gpmod
        def fset(self, value):
            """Set genomic prediction model"""
            check_is_GenomicModel(value, "gpmod")
            self._gpmod = value
        def fdel(self):
            """Delete genomic prediction model"""
            del self._gpmod
        return locals()
    gpmod = property(**gpmod())

    ################ Stochastic Parameters #################
    def var_err():
        doc = "Error variance for each trait."
        def fget(self):
            """Get error variance"""
            return numpy.repeat(1.0, self.gpmod.ntrait)
        def fset(self, value):
            """Set error variance"""
            error_readonly("var_err")
        def fdel(self):
            """Delete error variance"""
            error_readonly("var_err")
        return locals()
    var_err = property(**var_err())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def phenotype(self, pgmat, miscout = None, gpmod = None, **kwargs):
        """
        Phenotype a set of genotypes using a genomic prediction model.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes of the individuals to phenotype.
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        gpmod : GenomicModel, None
            Genomic prediction model to use to determine phenotypes.
            If None, use default genomic prediction model.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhenotypeDataFrame
            A PhenotypeDataFrame containing phenotypes for individuals.
        """
        # process arguments
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        if gpmod is None:
            gpmod = self.gpmod

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

    def set_h2(self, h2, pgmat, **kwargs):
        """
        Set the narrow sense heritability for environments.

        Parameters
        ----------
        h2 : float, numpy.ndarray
            Narrow sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise AttributeError("unsupported operation: heritability always set at 1.0")

    def set_H2(self, H2, pgmat, **kwargs):
        """
        Set the broad sense heritability for environments.

        Parameters
        ----------
        H2 : float, numpy.ndarray
            Broad sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise AttributeError("unsupported operation: heritability always set at 1.0")
