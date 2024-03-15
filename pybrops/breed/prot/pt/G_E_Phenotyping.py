"""
Module implementing phenotyping protocols for simulating phenotyping with no GxE
interaction.
"""

import copy
import math
from numbers import Integral, Real
from pathlib import Path
from typing import Optional, Union
import h5py
import numpy
from numpy.random import Generator, RandomState
import pandas

from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_type_python import check_is_Integral, check_is_dict
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group, check_h5py_File_is_readable, check_h5py_File_is_writable
from pybrops.core.error.error_value_numpy import check_ndarray_all_gt, check_ndarray_all_gteq, check_ndarray_ndim
from pybrops.core.random.prng import global_prng
from pybrops.breed.prot.pt.PhenotypingProtocol import PhenotypingProtocol
from pybrops.core.error.error_value_python import check_is_gt, check_is_gteq
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState, check_ndarray_dtype_is_integer, check_ndarray_dtype_is_real
from pybrops.core.error.error_value_numpy import check_ndarray_size
from pybrops.core.util.h5py import h5py_File_read_int, h5py_File_read_ndarray, h5py_File_read_ndarray_int, h5py_File_write_dict
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix

class G_E_Phenotyping(
        PhenotypingProtocol
    ):
    """
    Class implementing phenotyping protocols for simulating phenotyping with no
    GxE interaction.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            gpmod: GenomicModel, 
            nenv: Integral = 1, 
            nrep: Union[Integral,numpy.ndarray] = 1, 
            var_env: Optional[Union[Real,numpy.ndarray]] = None, 
            var_rep: Optional[Union[Real,numpy.ndarray]] = None, 
            var_err: Optional[Union[Real,numpy.ndarray]] = None, 
            rng: Optional[Union[Generator,RandomState]] = None, 
            **kwargs: dict
        ) -> None:
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

            If numeric, then broadcast ``var_err`` to an array of shape
            ``(ntrait,)``.

            If ``numpy.ndarray``, then must be of shape ``(ntrait,)``.
        
        rng : numpy.random.Generator, numpy.random.RandomState, None
            Random number source.

        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent initialization!
        self.gpmod = gpmod      # order 0
        self.nenv = nenv        # order 1
        self.nrep = nrep        # order 2
        self.var_env = var_env  # order 3
        self.var_rep = var_rep  # order 4
        self.var_err = var_err  # order 5
        self.rng = rng

    def __copy__(
            self
        ) -> 'G_E_Phenotyping':
        """
        Make a shallow copy of the ``G_E_Phenotyping`` object.

        Returns
        -------
        out : G_E_Phenotyping
            A shallow copy of the ``G_E_Phenotyping`` object.
        """
        # get the class
        cls = type(self)

        # create a new class object
        out = cls(
            gpmod   = copy.copy(self.gpmod),
            nenv    = copy.copy(self.nenv),
            nrep    = copy.copy(self.nrep),
            var_env = copy.copy(self.var_env),
            var_rep = copy.copy(self.var_rep),
            var_err = copy.copy(self.var_err),
            rng     = self.rng, # should not be copied
        )

        return out

    def __deepcopy__(
            self,
            memo: Optional[dict] = None,
        ) -> 'G_E_Phenotyping':
        """
        Make a deep copy of the ``G_E_Phenotyping`` object.

        Parameters
        ----------
        memo : dict, None
            An optional dictionary of memo metadata.

        Returns
        -------
        out : G_E_Phenotyping
            A deep copy of the ``G_E_Phenotyping`` object.
        """
        # get the class
        cls = type(self)

        # create a new class object
        out = cls(
            gpmod   = copy.deepcopy(self.gpmod, memo),
            nenv    = copy.deepcopy(self.nenv, memo),
            nrep    = copy.deepcopy(self.nrep, memo),
            var_env = copy.deepcopy(self.var_env, memo),
            var_rep = copy.deepcopy(self.var_rep, memo),
            var_err = copy.deepcopy(self.var_err, memo),
            rng     = self.rng, # should not be copied
        )

        return out

    ############################ Object Properties #############################

    ############### Genomic Model Properties ###############
    @property
    def gpmod(self) -> GenomicModel:
        """Genomic prediction model."""
        return self._gpmod
    @gpmod.setter
    def gpmod(self, value: GenomicModel) -> None:
        """Set genomic prediction model."""
        check_is_GenomicModel(value, "gpmod")
        self._gpmod = value

    ################ Replication Parameters ################
    @property
    def nenv(self) -> Integral:
        """Number of environments."""
        return self._nenv
    @nenv.setter
    def nenv(self, value: Integral) -> None:
        """Set number of environments"""
        check_is_Integral(value, "nenv")
        check_is_gt(value, "nenv", 0)
        self._nenv = value

    @property
    def nrep(self) -> numpy.ndarray:
        """Number of replications per environment."""
        return self._nrep
    @nrep.setter
    def nrep(self, value: Union[Integral,numpy.ndarray]) -> None:
        """Set number of replications per environment."""
        if isinstance(value, Integral):
            check_is_gt(value, "nrep", 0)
            value = numpy.full(self.nenv, value, int)
        elif isinstance(value, numpy.ndarray):
            check_ndarray_dtype_is_integer(value, "nrep")
            check_ndarray_ndim(value, "nrep", 1)
            check_ndarray_size(value, "nrep", self.nenv)
            check_ndarray_all_gt(value, "nrep", 0)
            self._nrep = value
        else:
            raise TypeError("'nrep' must be an integer type or numpy.ndarray")
        self._nrep = value

    ################ Stochastic Parameters #################
    @property
    def var_env(self) -> numpy.ndarray:
        """Variance across environments."""
        return self._var_env
    @var_env.setter
    def var_env(self, value: Union[Real,numpy.ndarray,None]) -> None:
        """Set variance across environments"""
        if value is None:
            value = numpy.full(self.gpmod.ntrait, 0.0, float)
        if isinstance(value, Real):
            check_is_gteq(value, "var_env", 0)                      # check var >= 0
            value = numpy.full(self.gpmod.ntrait, value, float)     # make env var matrix
        elif isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "var_env", 1)                 # check is 1d array
            check_ndarray_size(value, "var_env", self.gpmod.ntrait) # check sizes align
            check_ndarray_dtype_is_real(value, "var_env")           # check real dtype
            check_ndarray_all_gteq(value, "var_env", 0)             # check var >= 0
            if value.dtype != float:                                # if dtype != float64
                value = numpy.array(value, dtype = float)           # convert array to float64
        else:
            raise TypeError("'var_env' must be Real, numpy.ndarray, or None type")
        self._var_env = value

    @property
    def var_rep(self) -> numpy.ndarray:
        """Variance across replicates."""
        return self._var_rep
    @var_rep.setter
    def var_rep(self, value: Union[Real,numpy.ndarray]) -> None:
        """Set replicate variance"""
        if value is None:
            value = numpy.full(self.gpmod.ntrait, 0.0, float)
        if isinstance(value, Real):
            check_is_gteq(value, "var_rep", 0)                      # check var >= 0
            value = numpy.full(self.gpmod.ntrait, value, float)     # make rep var matrix
        elif isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "var_rep", 1)                 # check is 1d array
            check_ndarray_size(value, "var_rep", self.gpmod.ntrait) # check sizes align
            check_ndarray_dtype_is_real(value, "var_rep")           # check real dtype
            check_ndarray_all_gteq(value, "var_rep", 0)             # check var >= 0
            if value.dtype != float:                                # if dtype != float64
                value = numpy.array(value, dtype = float)           # convert array to float64
        else:
            raise TypeError("'var_rep' must be Real, numpy.ndarray, or None type")
        self._var_rep = value

    @property
    def var_err(self) -> numpy.ndarray:
        """Error variance for each trait."""
        return self._var_err
    @var_err.setter
    def var_err(self, value: Union[Real,numpy.ndarray]) -> None:
        """Set error variance"""
        if value is None:
            value = numpy.full(self.gpmod.ntrait, 0.0, float)
        if isinstance(value, Real):
            check_is_gteq(value, "var_err", 0)                      # check var >= 0
            value = numpy.full(self.gpmod.ntrait, value, float)     # make error var matrix
        elif isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "var_err", 1)                 # check is 1d array
            check_ndarray_size(value, "var_err", self.gpmod.ntrait) # check sizes align
            check_ndarray_dtype_is_real(value, "var_err")           # check real dtype
            check_ndarray_all_gteq(value, "var_err", 0)             # check var >= 0
            if value.dtype != float:                                # if dtype != float64
                value = numpy.array(value, dtype = float)           # convert array to float64
        else:
            raise TypeError("'var_err' must be Real, numpy.ndarray, or None type")
        self._var_err = value

    @property
    def rng(self) -> Union[Generator,RandomState]:
        """Random number generator."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[Generator,RandomState]) -> None:
        """Set random number generator."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng")
        self._rng = value

    ############################## Object Methods ##############################
    def copy(
            self
        ) -> 'G_E_Phenotyping':
        """
        Make a shallow copy of the ``G_E_Phenotyping`` object.

        Returns
        -------
        out : G_E_Phenotyping
            A shallow copy of the ``G_E_Phenotyping`` object.
        """
        return copy.copy(self)

    def deepcopy(
            self,
            memo: Optional[dict] = None,
        ) -> 'G_E_Phenotyping':
        """
        Make a deep copy of the ``G_E_Phenotyping`` object.

        Parameters
        ----------
        memo : dict, None
            An optional dictionary of memo metadata.

        Returns
        -------
        out : G_E_Phenotyping
            A deep copy of the ``G_E_Phenotyping`` object.
        """
        return copy.deepcopy(self, memo)
    
    def phenotype(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> pandas.DataFrame:
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : pandas.DataFrame
            A pandas.DataFrame containing phenotypes for individuals.
        """
        # check data types
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        if miscout is not None:
            check_is_dict(miscout, "miscout")

        # calculate true genotypic values
        gvmat = self.gpmod.gegv(pgmat)

        # (n,t) : calculate unscaled genotypic values
        mat = gvmat.unscale()

        # get the number of taxa in each replication
        ntaxa = gvmat.ntaxa

        # get the taxa names
        taxazfill = math.ceil(math.log10(gvmat.ntaxa))+1
        taxa = numpy.array(["Taxon"+str(i+1).zfill(taxazfill) for i in range(gvmat.ntaxa)], dtype=object) if gvmat.taxa is None else gvmat.taxa

        # get the taxa group names
        taxa_grp = None if gvmat.taxa_grp is None else gvmat.taxa_grp

        # lists to accumulate values
        taxa_ls = []
        taxa_grp_ls = None if gvmat.taxa_grp is None else []
        env_ls = []
        rep_ls = []
        values_ls = []

        # (t,) : environmental mean for traits
        env_mean = numpy.full(len(self.var_env), 0.0, float)

        # (t,t) : environmental covariance for traits
        env_cov = numpy.diag(self.var_env)

        # (t,) : replication mean for traits
        rep_mean = numpy.full(len(self.var_rep), 0.0, float)

        # (t,t) : replication covariance for traits
        rep_cov = numpy.diag(self.var_rep)

        # (t,) : error mean for traits
        err_mean = numpy.full(len(self.var_err), 0.0, float)

        # (t,t) : error covariance for traits
        err_cov = numpy.diag(self.var_err)

        # for each environment
        for env,env_nrep in zip(range(self.nenv),self.nrep):

            # (t,) : calculate the environmental effect
            env_effect = self.rng.multivariate_normal(env_mean, env_cov)

            # for each replicate in each environment
            for rep in range(env_nrep):

                # (t,) : calculate the replication effect
                rep_effect = self.rng.multivariate_normal(rep_mean, rep_cov)

                # (n,t) : calculate the error effect for each taxon
                err_effect = self.rng.multivariate_normal(err_mean, err_cov, ntaxa)

                # (n,t) : calculate phenotypic values for the rep within the env
                # (n,t) = (n,t) + (1,t) + (1,t) + (n,t)
                value = mat + env_effect[None,:] + rep_effect[None,:] + err_effect

                # append values
                taxa_ls.append(taxa)
                if taxa_grp_ls is not None:
                    taxa_grp_ls.append(taxa_grp)
                env_ls.append(numpy.repeat(env+1, ntaxa))
                rep_ls.append(numpy.repeat(rep+1, ntaxa))
                values_ls.append(value)
        
        # concatenate values
        # (nobs,)
        taxa_vt = numpy.concatenate(taxa_ls)
        # None or (nobs,)
        taxa_grp_vt = None if taxa_grp_ls is None else numpy.concatenate(taxa_grp_ls)
        # (nobs,)
        env_vt = numpy.concatenate(env_ls)
        # (nobs,)
        rep_vt = numpy.concatenate(rep_ls)
        # (nobs,t)
        values_mt = numpy.concatenate(values_ls, axis = 0)
        
        # construct dataframe labels
        labels_df = pandas.DataFrame({
            "taxa": taxa_vt,
            "taxa_grp": taxa_grp_vt,
            "env": env_vt,
            "rep": rep_vt
        })

        # construct dataframe data values
        traitzfill = math.ceil(math.log10(gvmat.ntrait))+1
        cols = ["Trait"+str(i+1).zfill(traitzfill) for i in range(gvmat.ntrait)] if gvmat.trait is None else gvmat.trait
        values_df = pandas.DataFrame(
            data = values_mt,
            columns = cols
        )

        # combine data labels and values
        out_df = pandas.concat([labels_df, values_df], axis = 1)

        return out_df

    def set_h2(
            self, 
            h2: Union[Real,numpy.ndarray], 
            pgmat: PhasedGenotypeMatrix, 
            **kwargs: dict
        ) -> None:
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
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")

        # get the additive variance of the founder population
        var_A = self.gpmod.var_A(pgmat)

        # calculate environmental variance
        # var_err = (1 - h2)/h2 * var_A - var_G
        # we assume var_G is zero, so var_err = (1 - h2)/h2 * var_A
        # scalar - (t,) -> (t,)
        # (t,) / (t,) -> (t,)
        # (t,) * (t,) -> (t,)
        self.var_err = (1.0 - h2) / h2 * var_A

    def set_H2(
            self, 
            H2: Union[Real,numpy.ndarray], 
            pgmat: PhasedGenotypeMatrix, 
            **kwargs: dict
        ) -> None:
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
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")

        # get the additive variance of the founder population
        var_G = self.gpmod.var_G(pgmat)

        # calculate environmental variance
        # var_err = (1 - h2)/h2 * var_A - var_G
        # we assume var_G is zero, so var_err = (1 - h2)/h2 * var_A
        # scalar - (t,) -> (t,)
        # (t,) / (t,) -> (t,)
        # (t,) * (t,) -> (t,)
        self.var_err = (1.0 - H2) / H2 * var_G

    ################### Matrix File I/O ####################
    def to_hdf5(
            self, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            overwrite: bool = True,
        ) -> None:
        """
        Write a ``G_E_Phenotyping`` object to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name to which to write. File is closed after writing.
            If ``h5py.File``, an already opened HDF5 file to which to write. File remains open after writing.

        groupname : str, None
            If ``str``, an HDF5 group name under which object data is stored.
            If ``None``, object is written to the base HDF5 group.

        overwrite : bool
            Whether to overwrite values in an HDF5 file if a field already exists.
        """
        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            # empty string
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        ############ process ``filename`` argument #############

        # h5 file object
        h5file = None
        
        # if we have a string or Path, open HDF5 file in append (``r+``) mode
        if isinstance(filename, (str,Path)):
            h5file = h5py.File(filename, "a")

        # elif we have an h5py.File, make sure mode is writable, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_writable(filename)
            h5file = filename

        # else raise TypeError
        else:
            raise TypeError(
                "filename must be of type ``str`` or ``h5py.File`` but received type {0}".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ### populate HDF5 file

        # data dictionary
        data = {
            "nenv"      : self.nenv,
            "nrep"      : self.nrep,
            "var_env"   : self.var_env,
            "var_rep"   : self.var_rep,
            "var_err"   : self.var_err,
        }

        # save data
        h5py_File_write_dict(h5file, groupname, data, overwrite)

        # close the file, only if the provided filename was a string or Path and not a h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            gpmod : GenomicModel = None,
        ) -> 'G_E_Phenotyping':
        """
        Read a ``G_E_Phenotyping`` from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which ``G_E_Phenotyping`` data is stored.
            If ``None``, ``G_E_Phenotyping`` is read from base HDF5 group.
        gpmod : GenomicModel
            A genomic model to bind to the ``G_E_Phenotyping`` protocol.
            This will be eliminated when a better storage mechanism is available.

        Returns
        -------
        out : G_E_Phenotyping
            An ``G_E_Phenotyping`` read from an HDF5 file.
        """
        ########################################################
        ############ process ``filename`` argument #############

        # HDF5 file object
        h5file = None

        # if we have a string or Path, open HDF5 file in append (``r``) mode
        if isinstance(filename, (str,Path)):
            check_file_exists(filename)
            h5file = h5py.File(filename, "r")

        # elif we have an h5py.File, make sure mode is in at least ``r`` mode, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_readable(filename)
            h5file = filename
        
        # else raise TypeError
        else:
            raise TypeError(
                "``filename`` must be of type ``str``, ``Path``, or ``h5py.File`` but received type ``{0}``".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # FIXME: errors if groupname == "" or "/"
            # if the group does not exist in the file, close and raise error
            check_h5py_File_has_group(h5file, groupname)

            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        ################## process ``gpmod`` ###################

        check_is_GenomicModel(gpmod, "gpmod")

        ########################################################
        ######## check that we have all required fields ########

        # all required arguments
        required_fields = []

        # for each required field, check if the field exists in the HDF5 file.
        for field in required_fields:
            check_h5py_File_has_group(h5file, groupname + field)
        
        ########################################################
        ### read data from HDF5 file and (optionally) close ####
        
        # output dictionary
        data = {
            "nenv"      : None,
            "nrep"      : None,
            "var_env"   : None,
            "var_rep"   : None,
            "var_err"   : None,
        }

        ##################################
        ### read mandatory data fields ###

        # read nenv data (dtype = int)
        data["nenv"] = h5py_File_read_int(h5file, groupname + "nenv")

        # read nrep array (ndarray dtype = int)
        data["nrep"] = h5py_File_read_ndarray_int(h5file, groupname + "nrep")

        #################################
        ### read optional data fields ###

        # read var_env array (ndarray dtype = any)
        if groupname + "var_env" in h5file:
            data["var_env"] = h5py_File_read_ndarray(h5file, groupname + "var_env")

        # read var_rep array (ndarray dtype = any)
        if groupname + "var_rep" in h5file:
            data["var_rep"] = h5py_File_read_ndarray(h5file, groupname + "var_rep")

        # read var_err array (ndarray dtype = any)
        if groupname + "var_err" in h5file:
            data["var_err"] = h5py_File_read_ndarray(h5file, groupname + "var_err")

        #####################################
        ### read optional metadata fields ###

        ######################
        ### close the file ###

        # close the file, only if the provided fieldname was a string or Path an not an h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

        ########################################################
        ### create object
        
        # create object from read data
        out = cls(
            gpmod   = gpmod,
            nenv    = data["nenv"],
            nrep    = data["nrep"],
            var_env = data["var_env"],
            var_rep = data["var_rep"],
            var_err = data["var_err"],
            rng     = None,
        )

        return out
