import numpy
import numbers
import pandas

import pybropt.core.random
from pybropt.breed.prot.pt.PhenotypingProtocol import PhenotypingProtocol
from pybropt.popgen.ptdf.DictPhenotypeDataFrame import DictPhenotypeDataFrame
from pybropt.model.gmod.GenomicModel import check_is_GenomicModel

from pybropt.core.error import check_is_positive
from pybropt.core.error import check_ndarray_is_1d
from pybropt.core.error import check_ndarray_size
from pybropt.core.error import check_is_Integral
from pybropt.core.error import check_ndarray_is_positive
from pybropt.core.error import check_ndarray_dtype_is_integer

class G_E_Phenotyping(PhenotypingProtocol):
    """docstring for G_E_Phenotyping."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, gpmod, nenv = 1, nrep = 1, var_env = 0, var_rep = 0, var_err = 0, rng = None, **kwargs):
        """
        Construct a phenotyping protocol that simulates environments as having
        a fixed effect, but no genotype by environment interaction. Variance
        across environments are equal.

        Parameters
        ----------
        gpmod : GenomicModel
            True genomic model to use for prediction.

        nenv : int
            Number of environments.
        nrep : int, numpy.ndarray
            Number of replications per environment.

            If ``int``, then broadcast ``nrep`` to an array of shape ``(nenv,)``
            If ``numpy.ndarray``, then must be of shape ``(nenv,)``
        var_env : numeric, numpy.ndarray
            Environmental variance parameter for each trait.
            Determines distribution of fixed effect added to each environment.

            If numeric, then broadcast ``var_env`` to an array of shape ``(ntrait,)``
            If ``numpy.ndarray``, then must be of shape ``(ntrait,)``
        var_rep : numeric, numpy.ndarray
            Replication variance parameter for each trait.

            Replication variance is assumed to be constant across environments.
            Replication is nested within each environment.

            If numeric, then broadcast ``var_rep`` to an array of shape ``(ntrait,)``
            If ``numpy.ndarray``, then must be of shape ``(ntrait,)``
        var_err : numeric, numpy.ndarray
            Error variance parameter.
            --------------------------------------------------------------------
            If numeric:
                Broadcast 'var_err' to an array of shape (ntrait,)
            If numpy.ndarray:
                Must be of shape (ntrait,)
            --------------------------------------------------------------------

            If numeric, then broadcast ``var_err`` to an array of shape ``(ntrait,)``
            If ``numpy.ndarray``, then must be of shape ``(ntrait,)``
        """
        super(G_E_Phenotyping, self).__init__(**kwargs)
        # order dependent initialization!
        self.gpmod = gpmod      # order 0
        self.nenv = nenv        # order 1
        self.nrep = nrep        # order 2
        self.var_env = var_env  # order 3
        self.var_rep = var_rep  # order 4
        self.var_err = var_err  # order 5
        self.rng = pybropt.core.random if rng is None else rng


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
    def var_env():
        doc = "Variance across environments"
        def fget(self):
            """Get variance across environments"""
            return self._var_env
        def fset(self, value):
            """Set variance across environments"""
            if isinstance(value, numbers.Number):
                check_is_positive(value, "var_env") # make sure >= 0 variance
                self._var_env = numpy.full(         # allocate empty array
                    self.gpmod.ntrait,              # ntrait length
                    value,                          # fill value
                    dtype = "float64"               # must be float64
                )
            elif isinstance(value, numpy.ndarray):
                check_ndarray_is_1d(value, "var_env")   # make sure is 1d array
                check_ndarray_size(                     # make sure size aligns with number of traits
                    value,
                    "var_env",
                    self.gpmod.ntrait
                )
                check_ndarray_is_positive(value, "var_env") # make sure we don't have negative variance
                if value.dtype != "float64":                # if dtype != float64
                    value = numpy.float64(value)            # convert array to float64
                self._var_env = value                       # set var_env
            else:
                raise TypeError("'var_env' must be a numeric or numpy.ndarray type")
        def fdel(self):
            """Delete variance across environments"""
            del self._var_env
        return locals()
    var_env = property(**var_env())

    def var_rep():
        doc = "Variance across replicates"
        def fget(self):
            """Get replicate variance"""
            return self._var_rep
        def fset(self, value):
            """Set replicate variance"""
            if isinstance(value, numbers.Number):
                check_is_positive(value, "var_rep") # make sure >= 0 variance
                self._var_rep = numpy.full(         # allocate empty array
                    self.gpmod.ntrait,              # ntrait length
                    value,                          # fill value
                    dtype = "float64"               # must be float64
                )
            elif isinstance(value, numpy.ndarray):
                check_ndarray_is_1d(value, "var_rep")   # make sure is 1d array
                check_ndarray_size(                     # make sure size aligns with number of traits
                    value,
                    "var_rep",
                    self.gpmod.ntrait
                )
                check_ndarray_is_positive(value, "var_rep") # make sure we don't have negative variance
                if value.dtype != "float64":                # if dtype != float64
                    value = numpy.float64(value)            # convert array to float64
                self._var_rep = value                       # set var_rep
            else:
                raise TypeError("'var_rep' must be a numeric or numpy.ndarray type")
            self._var_rep = value
        def fdel(self):
            """Delete replicate variance"""
            del self._var_rep
        return locals()
    var_rep = property(**var_rep())

    def var_err():
        doc = "Error variance for each trait."
        def fget(self):
            """Get error variance"""
            return self._var_err
        def fset(self, value):
            """Set error variance"""
            if isinstance(value, numbers.Number):
                check_is_positive(value, "var_err") # make sure >= 0 variance
                self._var_err = numpy.full(         # allocate empty array
                    self.gpmod.ntrait,              # ntrait length
                    value,                          # fill value
                    dtype = "float64"               # must be float64
                )
            elif isinstance(value, numpy.ndarray):
                check_ndarray_is_1d(value, "var_err")   # make sure is 1d array
                check_ndarray_size(                     # make sure size aligns with number of traits
                    value,
                    "var_err",
                    self.gpmod.ntrait
                )
                check_ndarray_is_positive(value, "var_err") # make sure we don't have negative variance
                if value.dtype != "float64":                # if dtype != float64
                    value = numpy.float64(value)            # convert array to float64
                self._var_err = value                       # set var_err
            else:
                raise TypeError("'var_err' must be a numeric or numpy.ndarray type")
        def fdel(self):
            """Delete error variance"""
            del self._var_err
        return locals()
    var_err = property(**var_err())

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
    def phenotype(self, pgmat, miscout = None, gpmod = None, nenv = None, nrep = None, var_env = None, var_rep = None, var_err = None, **kwargs):
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
        nenv : int
            Number of environments.
        nrep : int, numpy.ndarray
            Number of replications per environment.
        var_env : numeric, numpy.ndarray
            Environmental variance parameter.
        var_rep : numeric, numpy.ndarray
            Replication variance parameter.
        var_err : numeric, numpy.ndarray
            Error variance parameter.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhenotypeDataFrame
            A PhenotypeDataFrame containing phenotypes for individuals.
        """
        ################### set default parameters if needed ###################
        # set default gpmod
        if gpmod is None:
            gpmod = self.gpmod
        else:
            check_is_GenomicModel(gpmod, "gpmod")

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

        # set default var_env
        if var_env is None:
            var_env = self.var_env
        elif isinstance(var_env, numbers.Number):   # if var_env is a number
            check_is_positive(var_env, "var_env")   # make sure >= 0 variance
            var_env = numpy.full(                   # allocate empty array
                gpmod.ntrait,                       # ntrait length
                var_env,                            # fill value
                dtype = "float64"                   # must be float64
            )
        elif isinstance(var_env, numpy.ndarray):                # if var_env is a numpy.ndarray
            check_ndarray_is_1d(var_env, "var_env")             # make sure is 1d array
            check_ndarray_size(var_env, "var_env", gpmod.ntrait)# make sure size aligns with number of traits
            check_ndarray_is_positive(var_env, "var_env")       # make sure we don't have negative variance
            if var_env.dtype != "float64":                      # if dtype != float64
                var_env = numpy.float64(var_env)                # convert array to float64
        else:
            raise TypeError("'var_env' must be a numeric or numpy.ndarray type")

        # set default var_rep
        if var_rep is None:
            var_rep = self.var_rep
        elif isinstance(var_rep, numbers.Number):   # if var_rep is a number
            check_is_positive(var_rep, "var_rep")   # make sure >= 0 variance
            var_rep = numpy.full(                   # allocate empty array
                gpmod.ntrait,                       # ntrait length
                var_rep,                            # fill value
                dtype = "float64"                   # must be float64
            )
        elif isinstance(var_rep, numpy.ndarray):                # if var_rep is a numpy.ndarray
            check_ndarray_is_1d(var_rep, "var_rep")             # make sure is 1d array
            check_ndarray_size(var_rep, "var_rep", gpmod.ntrait)# make sure size aligns with number of traits
            check_ndarray_is_positive(var_rep, "var_rep")       # make sure we don't have negative variance
            if var_rep.dtype != "float64":                      # if dtype != float64
                var_rep = numpy.float64(var_rep)                # convert array to float64
        else:
            raise TypeError("'var_rep' must be a numeric or numpy.ndarray type")

        # set default var_err
        if var_err is None:
            var_err = self.var_err
        elif isinstance(var_err, numbers.Number):   # if var_err is a number
            check_is_positive(var_err, "var_err")   # make sure >= 0 variance
            var_err = numpy.full(                   # allocate empty array
                gpmod.ntrait,                       # ntrait length
                var_err,                            # fill value
                dtype = "float64"                   # must be float64
            )
        elif isinstance(var_err, numpy.ndarray):                # if var_err is a numpy.ndarray
            check_ndarray_is_1d(var_err, "var_err")             # make sure is 1d array
            check_ndarray_size(var_err, "var_err", gpmod.ntrait)# make sure size aligns with number of traits
            check_ndarray_is_positive(var_err, "var_err")       # make sure we don't have negative variance
            if var_err.dtype != "float64":                      # if dtype != float64
                var_err = numpy.float64(var_err)                # convert array to float64
        else:
            raise TypeError("'var_err' must be a numeric or numpy.ndarray type")

        ################## extract breeding value information ##################
        bvmat = gpmod.gebv(pgmat)   # gather true breeding values
        mat = bvmat.descale()       # get descaled breeding values
        taxa = bvmat.taxa           # taxa names
        taxa_grp = bvmat.taxa_grp   # taxa groups
        trait = bvmat.trait         # trait names
        ntaxa = bvmat.ntaxa         # number of taxa
        ntrait = bvmat.ntrait       # number of traits
        nplot = nrep.sum()          # get number of plots per genotype

        ############### perform error checks ###############
        # make sure all taxa names are all valid
        if taxa is None:
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' does not have taxa names")
        if any(e is None for e in taxa):
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' has taxa name(s) which are 'None'")

        # make sure all trait names are valid
        if trait is None:
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' does not have trait names")
        if any(e is None for e in trait):
            raise ValueError("unable to construct phenotype dataframe: breeding value matrix produced by 'gpmod.predict(pgmat)' has trait name(s) which are 'None'")

        ############################ construct data ############################
        # allocate memory for phenotypes
        pheno_mat = numpy.empty((ntaxa * nplot, ntrait), dtype = "float64")

        # copy breeding value information into pheno_mat
        for i in range(nplot):
            pheno_mat[i*ntaxa:(i+1)*ntaxa,:] = mat

        # sample environmental effects (additive; no GxE)
        env_effect = self.rng.normal(
            0.0,                    # sample mean = 0.0
            numpy.sqrt(var_env),    # standard deviation
            (nenv, ntrait)          # sample nenv times
        )

        # add environmental effects to the breeding values
        stix = ntaxa * (nrep.cumsum() - nrep)       # calculate start indices
        spix = ntaxa * (nrep.cumsum())              # calculate stop indices
        for i,(st,sp) in enumerate(zip(stix,spix)): # for each environment block
            pheno_mat[st:sp,:] += env_effect[i,:]   # add environmental effect

        # sample replicate effects (additive; no interaction)
        rep_effect = self.rng.normal(
            0.0,                    # sample mean = 0.0
            numpy.sqrt(var_rep),    # standard deviation
            (nplot, ntrait)         # sample nplot times
        )

        # add replicate effects to the breeding values
        stix = range(0, ntaxa * nplot, ntaxa)           # calculate start indices
        spix = range(ntaxa, ntaxa * (nplot + 1), ntaxa) # calculate stop indices
        for i,(st,sp) in enumerate(zip(stix,spix)):     # for each environment block
            pheno_mat[st:sp,:] += rep_effect[i,:]       # add environmental effect

        # sample error effects
        err_effect = self.rng.normal(
            0.0,                    # sample mean = 0.0
            numpy.sqrt(var_err),    # sample standard deviation
            pheno_mat.shape         # (ntaxa * nplot, ntrait)
        )

        # add error effects to the breeding values
        pheno_mat += err_effect

        # create taxa vector
        taxa_vec = numpy.tile(taxa, nplot)

        # create taxa_grp vector
        taxa_grp_vec = None if taxa_grp is None else numpy.tile(taxa_grp, nplot)

        # create environment vector
        env_vec = numpy.repeat(numpy.repeat(numpy.arange(nenv), nrep), ntaxa)

        # create replicate vector
        rep_vec = numpy.repeat(numpy.concatenate([numpy.arange(e) for e in nrep]), ntaxa)

        ################ construct data dictionary and metadata ################
        data = {}                   # dictionary to contain all data
        col_analysis_type = {}
        col_analysis_effect = {}
        analysis_type = []      # list to contain variable type metadata
        analysis_effect = []    # list to contain variable effect type metadata


        numpy_to_R_dtype_lookup = {
            numpy.dtype('float16'): "double",
            numpy.dtype('float32'): "double",
            numpy.dtype('float64'): "double",
            numpy.dtype('int8'): "int",
            numpy.dtype('int16'): "int",
            numpy.dtype('int32'): "int",
            numpy.dtype('int64'): "int",
            numpy.dtype('uint8'): "int",
            numpy.dtype('uint16'): "int",
            numpy.dtype('uint32'): "int",
            numpy.dtype('uint64'): "int",
            numpy.dtype('object'): "factor(str)",
            numpy.dtype('bool'): "bool"
        }

        for i,t in enumerate(trait):
            colview = pheno_mat[:,i]
            data[t] = colview
            try:
                col_analysis_type[t] = numpy_to_R_dtype_lookup[colview.dtype]
            except KeyError:
                raise TypeError("unrecognized type {0}".format(colview.dtype))
            col_analysis_effect[t] = "response"

        # add taxa fixed effect
        data["taxa"] = taxa_vec                     # add taxa to dict
        col_analysis_type["taxa"] = "factor(str)"   # taxa is a string factor
        col_analysis_effect["taxa"] = "fixed"       # taxa is a fixed effect

        # add taxa group (family) fixed effect
        data["taxa_grp"] = taxa_grp_vec
        col_analysis_type["taxa_grp"] = None if taxa_grp_vec is None else "factor(int)"
        col_analysis_effect["taxa_grp"] = None if taxa_grp_vec is None else "fixed"

        # add environment random effect
        data["env"] = env_vec                   # add environment to dict
        col_analysis_type["env"] = "factor(int)"# environment is an int factor
        col_analysis_effect["env"] = "random"   # environment is a random effect

        # add replicate random effect
        data["rep"] = rep_vec                       # add replication to dict
        col_analysis_type["rep"] = "factor(int)"    # replication is an int factor
        col_analysis_effect["rep"] = "random"       # replication is a random effect

        # construct DictPhenotypeDataFrame
        ptdf = DictPhenotypeDataFrame(
            data = data,
            col_analysis_type = col_analysis_type,
            col_analysis_effect = col_analysis_effect,
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
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        # set default gpmod
        if gpmod is None:
            gpmod = self.gpmod
        else:
            check_is_GenomicModel(gpmod, "gpmod")

        # get GEBVs
        gebv = gpmod.gebv(pgmat)

        # get variance of breeding values
        var_A = gebv.tvar(descale = True) # (t,)

        # calculate environmental variance
        # var_err = (1 - h2)/h2 * var_A - var_G
        # we assume var_G is zero, so var_err = (1 - h2)/h2 * var_A
        # scalar - (t,) -> (t,)
        # (t,) / (t,) -> (t,)
        # (t,) * (t,) -> (t,)
        self.var_err = (1.0 - h2) / h2 * var_A

    def set_H2(self, H2, pgmat, gpmod = None, **kwargs):
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
        raise NotImplementedError("method is abstract")
