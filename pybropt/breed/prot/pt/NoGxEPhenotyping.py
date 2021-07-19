import numpy
import numbers
import pandas

from . import PhenotypingProtocol
from pybropt.popgen.ptdf import PandasPhenotypeDataFrame
from pybropt.model.gmod import check_is_GenomicModel

from pybropt.core.error import check_is_positive
from pybropt.core.error import check_ndarray_is_1d
from pybropt.core.error import check_ndarray_size


class NoGxEPhenotyping(PhenotypingProtocol):
    """docstring for NoGxEPhenotyping."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, gpmod, var_E = 0, nenv = 1, nrep = 1, **kwargs):
        """
        Parameters
        ----------
        gpmod : GenomicModel
            True genomic model to use for prediction.
        var_E : numeric, numpy.ndarray
            Environmental variance parameter.
        nenv : int
            Number of environments.
        nrep : int, numpy.ndarray
            Number of replications per environment.
        """
        super(NoGxEPhenotyping, self).__init__(**kwargs)
        # order dependent initialization!
        self.gpmod = gpmod  # order 0
        self.var_E = var_E  # order 1
        self.nenv = nenv    # order 2
        self.nrep = nrep    # order 3

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
    def var_E():
        doc = "Environmental variance for each trait."
        def fget(self):
            """Get environmental variance"""
            return self._var_E
        def fset(self, value):
            """Set environmental variance"""
            ntrait = self.gpmod.ntrait
            if isinstance(value, numbers.Number):
                check_is_positive(value, "var_E")   # make sure >= 0 variance
                self._var_E = numpy.full(           # allocate empty array
                    ntrait,                         # ntrait length
                    var_E,                          # fill value
                    dtype = "float64"               # must be float64
                )
            elif isinstance(value, numpy.ndarray):
                check_ndarray_is_1d(value, "var_E", 1)      # make sure is 1d array
                check_ndarray_size(value, "var_E", ntrait)  # make sure size aligns with number of traits
                check_ndarray_is_positive(value, "var_E")   # make sure we don't have negative variance
                if value.dtype != "float64":                # if dtype != float64
                    value = numpy.float64(value)            # convert array to float64
                self._var_E = value                         # set var_E
            else:
                raise TypeError("'var_E' must be a numeric or numpy.ndarray type")
        def fdel(self):
            """Delete environmental variance"""
            del self._var_E
        return locals()
    var_E = property(**var_E())

    ################ Replication Parameters ################
    def nenv():
        doc = "Number of environments"
        def fget(self):
            """Get number of environments"""
            return self._nenv
        def fset(self, value):
            """Set number of environments"""
            check_is_Integral(value, "nenv")
            self._nenv = value
        def fdel(self):
            """Delete number of environments"""
            del self._nenv
        return locals()
    nenv = property(**nenv())

    def nrep():
        doc = "Number of replications per environment"
        def fget(self):
            """Get number of replications per environment"""
            return self._nrep
        def fset(self, value):
            """Set number of replications per environment"""
            if isinstance(value, numbers.Integral):
                check_is_positive(value, "nrep")
                self._nrep = numpy.repeat(value, self.nenv)
            elif isinstance(value, numpy.ndarray):
                check_ndarray_dtype_is_integer(value, "nrep")
                check_ndarray_is_1d(value, "nrep")
                check_ndarray_is_positive(value, "nrep")
                check_ndarray_size(value, "nrep", self.nenv)
                self._nrep = value
            else:
                raise TypeError("'nrep' must be an integer type or numpy.ndarray")
        def fdel(self):
            """Delete number of replications per environment"""
            del self._nrep
        return locals()
    nrep = property(**nrep())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def phenotype(self, pgmat, gpmod = None, var_E = None, nenv = None, nrep = None, **kwargs):
        """
        Phenotype a set of genotypes using a genomic prediction model.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes of the individuals to phenotype.
        gpmod : GenomicModel, None
            Genomic prediction model to use to determine phenotypes.
        var_E : numeric, numpy.ndarray
            Environmental variance parameter.
        nenv : int
            Number of environments.
        nrep : int, numpy.ndarray
            Number of replications per environment.
        **kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhenotypeDataFrame
            DataFrame containing phenotypes.
        """
        ################### set default parameters if needed ###################
        # set default gpmod
        if gpmod is None:
            gpmod = self.gpmod
        else:
            check_is_GenomicModel(gpmod, "gpmod")

        # set default var_E
        if var_E is None:
            var_E = self.var_E
        elif isinstance(var_E, numbers.Number): # if var_E is a number
            check_is_positive(var_E, "var_E")   # make sure >= 0 variance
            var_E = numpy.full(                 # allocate empty array
                gpmod.ntrait,                   # ntrait length
                var_E,                          # fill value
                dtype = "float64"               # must be float64
            )
        elif isinstance(var_E, numpy.ndarray):                  # if var_E is a numpy.ndarray
            check_ndarray_is_1d(var_E, "var_E")                 # make sure is 1d array
            check_ndarray_size(var_E, "var_E", gpmod.ntrait)    # make sure size aligns with number of traits
            check_ndarray_is_positive(var_E, "var_E")           # make sure we don't have negative variance
            if var_E.dtype != "float64":                        # if dtype != float64
                var_E = numpy.float64(var_E)                    # convert array to float64
        else:
            raise TypeError("'var_E' must be a numeric or numpy.ndarray type")

        # set default nenv
        if nenv is None:
            nenv = self.nenv                    # set to default nenv
        else:
            check_is_Integral(nenv, "nenv")     # error check nenv

        # set default nrep
        if nrep is None:
            nrep = self.nrep
        elif isinstance(nrep, numbers.Integral):
            check_is_positive(nrep, "nrep")
            nrep = numpy.repeat(nrep, nenv)
        elif isinstance(nrep, numpy.ndarray):
            check_ndarray_dtype_is_integer(nrep, "nrep")
            check_ndarray_is_1d(nrep, "nrep")
            check_ndarray_is_positive(nrep, "nrep")
            check_ndarray_size(nrep, "nrep", nenv)
            nrep = value
        else:
            raise TypeError("'nrep' must be an integer type or numpy.ndarray")

        ################## extract breeding value information ##################
        # gather true breeding values
        bvmat = gpmod.predict(pgmat)

        # gather pointers to raw matrices and parameters
        mat = bvmat.mat             # breeding values
        taxa = bvmat.taxa           # taxa names
        taxa_grp = bvmat.taxa_grp   # taxa groups
        trait = bvmat.trait         # trait names
        ndim = bvmat.mat_ndim       # number of matrix dimensions
        taxis = bvmat.trait_axis    # trait axis
        nplot = nrep.sum()          # get number of plots per genotype

        # perform error checks
        if taxa is None:
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' does not have taxa names")
        if any(e is None for e in taxa):
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' has taxa name(s) which are 'None'")
        if trait is None:
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' does not have trait names")
        if any(e is None for e in trait):
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' has trait name(s) which are 'None'")
        if ndim != 2:
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix is not 2d.")

        ############################ construct data ############################
        # construct dictionary
        data_dict = {}

        # add taxa names
        data_dict["taxa"] = numpy.tile(taxa, nplot) # tile genotypes for number of plots per genotype

        # if there are taxa groups, add group information
        if taxa_grp is not None:
            data_dict["taxa_grp"] = numpy.tile(taxa_grp, nplot)

        # add each trait and corresponding data
        for i,(e,t) in enumerate(zip(var_E, trait)):
            s = tuple(i if a == taxis else slice(None) for a in range(ndim))    # construct matrix slice selection tuple
            traitpt = numpy.tile(mat[s], nplot)                                 # tile true breeding values
            traitpt = traitpt + self.rng.normal(0.0,numpy.sqrt(e),len(traitpt)) # add noise to phenotypes
            data_dict[t] = traitpt                                              # select data and put into dictionary

        # construct pandas.DataFrame
        df = pandas.DataFrame(data_dict)

        # construct PandasPhenotypeDataFrame
        ptdf = PandasPhenotypeDataFrame(
            df = df,
            **kwargs
        )

        return ptdf

    def set_h2(self, h2, pgmat, gpmod = None, **kwargs):
        """
        Set the narrow sense heritability for environments.

        Parameters
        ----------
        h2 : float, numpy.ndarray
            Narrow sense heritability.
        pgmat : PhasedGenotypeVariantMatrix
            Founder genotypes.
        **kwargs : dict
            Additional keyword arguments
        """
        # set default gpmod
        if gpmod is None:
            gpmod = self.gpmod
        else:
            check_is_GenomicModel(gpmod, "gpmod")

        # get matrices
        x = pgmat.tacount() # (n,p) genotype matrix
        b = gpmod.beta      # (p,t) marker effect matrix

        # get breeding values
        y = x @ b           # (n,p) @ (p,t) -> (n,t)

        # get variance of breeding values
        var_A = y.var(0)    # (n,t) -> (t,)

        # calculate environmental variance
        # var_E = (1 - h2)/h2 * var_A - var_G
        # we assume var_G is zero, so var_E = (1 - h2)/h2 * var_A
        # scalar - (t,) -> (t,)
        # (t,) / (t,) -> (t,)
        # (t,) * (t,) -> (t,)
        self.var_E = (1.0 - h2) / h2 * var_A

    def set_H2(self, H2, pgmat, gpmod = None, **kwargs):
        """
        Set the broad sense heritability for environments.

        Parameters
        ----------
        H2 : float, numpy.ndarray
            Broad sense heritability.
        pgmat : PhasedGenotypeVariantMatrix
            Founder genotypes.
        **kwargs : dict
            Additional keyword arguments
        """
        raise NotImplementedError("method is abstract")
