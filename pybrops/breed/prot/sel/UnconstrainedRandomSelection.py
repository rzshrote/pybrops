"""
Module implementing selection protocols for random selection.
"""

from typing import Callable, Union
import numpy
import types
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.breed.prot.sel.UnconstrainedSelectionProtocol import UnconstrainedSelectionProtocol
from pybrops.core.error import check_is_bool
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_callable
from pybrops.core.random.prng import global_prng

class RandomSelection(UnconstrainedSelectionProtocol):
    """
    Class implementing selection protocols for random selection.

    # TODO: add formulae for methodology.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            nparent: int, 
            ncross: int, 
            nprogeny: int, 
            replace: bool = False, 
            objfn_trans = None, 
            objfn_trans_kwargs = None, 
            objfn_wt = 1.0, 
            ndset_trans = None, 
            ndset_trans_kwargs = None, 
            ndset_wt = 1.0, 
            rng = None, 
            **kwargs: dict
        ):
        """
        Construct a RandomSelection operator.

        Parameters
        ----------
        nparent : int
            Number of parents to select.
        ncross : int
            Number of crosses per configuration.
        nprogeny : int
            Number of progeny to derive from each cross.
        replace : bool
            Whether to sample parents with replacement.
        rng : numpy.random.Generator or None
            A random number generator source. Used for optimization algorithms.
            If ``rng`` is ``None``, use ``pybrops.core.random`` module
            (NOT THREAD SAFE!).
        """
        super(RandomSelection, self).__init__(**kwargs)

        # assign variables
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.replace = replace
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = objfn_trans_kwargs
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs
        self.ndset_wt = ndset_wt
        self.rng = rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def nparent(self) -> int:
        """Number of parents to select."""
        return self._nparent
    @nparent.setter
    def nparent(self, value: int) -> None:
        """Set number of parents to select."""
        check_is_int(value, "nparent")      # must be int
        check_is_gt(value, "nparent", 0)    # int must be >0
        self._nparent = value
    @nparent.deleter
    def nparent(self) -> None:
        """Delete number of parents to select."""
        del self._nparent

    @property
    def ncross(self) -> int:
        """Number of crosses per configuration."""
        return self._ncross
    @ncross.setter
    def ncross(self, value: int) -> None:
        """Set number of crosses per configuration."""
        check_is_int(value, "ncross")       # must be int
        check_is_gt(value, "ncross", 0)     # int must be >0
        self._ncross = value
    @ncross.deleter
    def ncross(self) -> None:
        """Delete number of crosses per configuration."""
        del self._ncross

    @property
    def nprogeny(self) -> int:
        """Number of progeny to derive from each cross configuration."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: int) -> None:
        """Set number of progeny to derive from each cross configuration."""
        check_is_int(value, "nprogeny")     # must be int
        check_is_gt(value, "nprogeny", 0)   # int must be >0
        self._nprogeny = value
    @nprogeny.deleter
    def nprogeny(self) -> None:
        """Delete number of progeny to derive from each cross configuration."""
        del self._nprogeny

    @property
    def replace(self) -> bool:
        """Whether to replace individuals when sampling."""
        return self._replace
    @replace.setter
    def replace(self, value: bool) -> None:
        """Set whether to replace individuals when sampling."""
        check_is_bool(value, "replace")
        self._replace = value
    @replace.deleter
    def replace(self) -> None:
        """Delete whether to replace individuals when sampling."""
        del self._replace

    @property
    def objfn_trans(self) -> Union[Callable,None]:
        """Objective function transformation function."""
        return self._objfn_trans
    @objfn_trans.setter
    def objfn_trans(self, value: Union[Callable,None]) -> None:
        """Set objective function transformation function."""
        if value is not None:                       # if given object
            check_is_callable(value, "objfn_trans") # must be callable
        self._objfn_trans = value
    @objfn_trans.deleter
    def objfn_trans(self) -> None:
        """Delete objective function transformation function."""
        del self._objfn_trans

    @property
    def objfn_trans_kwargs(self) -> dict:
        """Objective function transformation function keyword arguments."""
        return self._objfn_trans_kwargs
    @objfn_trans_kwargs.setter
    def objfn_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set objective function transformation function keyword arguments."""
        if value is None:                           # if given None
            value = {}                              # set default to empty dict
        check_is_dict(value, "objfn_trans_kwargs")  # check is dict
        self._objfn_trans_kwargs = value
    @objfn_trans_kwargs.deleter
    def objfn_trans_kwargs(self) -> None:
        """Delete objective function transformation function keyword arguments."""
        del self._objfn_trans_kwargs

    # TODO: finish error checks
    @property
    def objfn_wt(self) -> Union[float,numpy.ndarray]:
        """Objective function weights."""
        return self._objfn_wt
    @objfn_wt.setter
    def objfn_wt(self, value: Union[float,numpy.ndarray]) -> None:
        """Set objective function weights."""
        self._objfn_wt = value
    @objfn_wt.deleter
    def objfn_wt(self) -> None:
        """Delete objective function weights."""
        del self._objfn_wt

    @property
    def ndset_trans(self) -> Union[Callable,None]:
        """Nondominated set transformation function."""
        return self._ndset_trans
    @ndset_trans.setter
    def ndset_trans(self, value: Union[Callable,None]) -> None:
        """Set nondominated set transformation function."""
        if value is not None:                       # if given object
            check_is_callable(value, "ndset_trans") # must be callable
        self._ndset_trans = value
    @ndset_trans.deleter
    def ndset_trans(self) -> None:
        """Delete nondominated set transformation function."""
        del self._ndset_trans

    @property
    def ndset_trans_kwargs(self) -> dict:
        """Nondominated set transformation function keyword arguments."""
        return self._ndset_trans_kwargs
    @ndset_trans_kwargs.setter
    def ndset_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set nondominated set transformation function keyword arguments."""
        if value is None:                           # if given None
            value = {}                              # set default to empty dict
        check_is_dict(value, "ndset_trans_kwargs")  # check is dict
        self._ndset_trans_kwargs = value
    @ndset_trans_kwargs.deleter
    def ndset_trans_kwargs(self) -> None:
        """Delete nondominated set transformation function keyword arguments."""
        del self._ndset_trans_kwargs

    # TODO: finish error checks
    @property
    def ndset_wt(self) -> Union[float,numpy.ndarray]:
        """Nondominated set weights."""
        return self._ndset_wt
    @ndset_wt.setter
    def ndset_wt(self, value: Union[float,numpy.ndarray]) -> None:
        """Set nondominated set weights."""
        self._ndset_wt = value
    @ndset_wt.deleter
    def ndset_wt(self) -> None:
        """Delete nondominated set weights."""
        del self._ndset_wt

    @property
    def rng(self) -> Union[numpy.random.Generator,numpy.random.RandomState]:
        """Random number generator source."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[numpy.random.Generator,numpy.random.RandomState]) -> None:
        """Set random number generator source."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng") # check is numpy.Generator
        self._rng = value
    @rng.deleter
    def rng(self) -> None:
        """Delete random number generator source."""
        del self._rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, nparent = None, ncross = None, nprogeny = None, replace = None, **kwargs: dict):
        """
        Select parents for breeding.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Phased genotype matrix containing full genome information.
        gmat : GenotypeMatrix
            Genotype matrix containing genotype data (phased or unphased)
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        method : str
            Options: "single", "pareto"
        nparent : int, None
            Number of parents. If ``None``, use default.
        ncross : int, None
            Number of crosses per configuration. If ``None``, use default.
        nprogeny : int
            Number of progeny per cross. If ``None``, use default.
        replace : bool, None
            Whether to sample parents with or without replacement. If ``None``,
            use default.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing four objects: ``(pgmat, sel, ncross, nprogeny)``.

            Where:

            - ``pgmat`` is a PhasedGenotypeMatrix of parental candidates.
            - ``sel`` is a ``numpy.ndarray`` of indices specifying a cross
              pattern. Each index corresponds to an individual in ``pgmat``.
            - ``ncross`` is a ``numpy.ndarray`` specifying the number of
              crosses to perform per cross pattern.
            - ``nprogeny`` is a ``numpy.ndarray`` specifying the number of
              progeny to generate per cross.
        """
        # get default parameters if any are None
        if nparent is None:
            nparent = self.nparent
        if ncross is None:
            ncross = self.ncross
        if nprogeny is None:
            nprogeny = self.nprogeny
        if replace is None:
            replace = self.replace

        # make random selection decision
        sel = self.rng.choice(
            pgmat.ntaxa,            # number of taxa to select from
            size = nparent,
            replace = replace
        )

        return pgmat, sel, ncross, nprogeny

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, **kwargs: dict):
        """
        Return a selection objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Not used by this function.
        gmat : GenotypeMatrix
            Input genotype matrix.
        ptdf : PhenotypeDataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Not used by this function.
        gpmod : LinearGenomicModel
            Linear genomic prediction model.

        Returns
        -------
        outfn : function
            A selection objective function for the specified problem.
        """
        # get default parameters if any are None
        if trans is None:
            trans = self.objfn_trans
        if trans_kwargs is None:
            trans_kwargs = self.objfn_trans_kwargs

        # get parameters/pointers
        n = gmat.ntaxa      # get number of taxa
        t = gpmod.ntrait    # get number of traits
        rng = self.rng      # get random number generator

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,         # byte code pointer
            self.objfn_static.__globals__,      # global variables
            None,                               # new name for the function
            (n, t, rng, trans, trans_kwargs),   # default values for arguments
            self.objfn_static.__closure__       # closure byte code pointer
        )

        return outfn

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, **kwargs: dict):
        """
        Return a vectorized selection objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Not used by this function.
        gmat : GenotypeMatrix
            Input genotype matrix.
        ptdf : PhenotypeDataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Not used by this function.
        gpmod : LinearGenomicModel
            Linear genomic prediction model.

        Returns
        -------
        outfn : function
            A vectorized selection objective function for the specified problem.
        """
        # get default parameters if any are None
        if trans is None:
            trans = self.objfn_trans
        if trans_kwargs is None:
            trans_kwargs = self.objfn_trans_kwargs

        # get parameters/pointers
        n = gmat.mat.shape[0]   # get number of taxa
        t = gpmod.beta.shape[1] # get number of traits
        rng = self.rng          # get random number generator

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,     # byte code pointer
            self.objfn_vec_static.__globals__,  # global variables
            None,                               # new name for the function
            (n, t, rng, trans, trans_kwargs),   # default values for arguments
            self.objfn_vec_static.__closure__   # closure byte code pointer
        )

        return outfn

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs: dict):
        """
        Random selection has no Pareto frontier since it has no objective function.
        Raises RuntimeError.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.
        """
        raise RuntimeError("Random selection has no Pareto frontier since it has no objective function.")

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, n, t, rng, trans, kwargs):
        """
        Randomly assign a score in the range :math:`[-1,1)` individuals for each trait.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            If ``sel`` is None, use all individuals.
        n : int
            The number of individuals available for selection.
        t : int
            The number of traits.
        rng : numpy.random.Generator, numpy.random.RandomState
            Random number generator.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        cgs : numpy.ndarray
            A GEBV matrix of shape ``(t,)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.
        """
        # get random number vector lengths
        nsel = n if sel is None else len(sel)   # get number of individuals
        ntrait = t                              # get number of traits

        # randomly sample from uniform distribution
        # Step 1: (k,t)                 # randomly score individuals
        # Step 2: (k,t).sum(0) -> (t,)  # sum across individuals
        rs = rng.uniform(-1.0, 1.0, (nsel,ntrait)).sum(0)

        # apply transformations
        # (t,) ---trans---> (?,)
        if trans:
            rs = trans(rs, **kwargs)

        return rs

    @staticmethod
    def objfn_vec_static(sel, n, t, rng, trans, kwargs):
        """
        Randomly assign a score in the range :math:`[-1,1)` individuals for
        each trait.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape ``(j,k)``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            If ``sel`` is None, score each individual separately: (n,1)
        n : int
            The number of individuals available for selection.
        t : int
            The number of traits.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single numpy.ndarray argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        cgs : numpy.ndarray
            A GEBV matrix of shape ``(j,t)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``t`` is the number of traits.
        """
        # get random number vector lengths
        nsel = (n,1) if sel is None else sel.shape  # get number of individuals (j,k)
        ntrait = (t,)                               # get number of traits (t,)

        # randomly sample from uniform distribution
        # Step 1: (j,k,t)                   # randomly score individuals
        # Step 2: (j,k,t).sum(1) -> (j,t)   # sum across individuals
        rs = rng.uniform(-1.0, 1.0, nsel+ntrait).sum(1)

        # apply transformations
        # (j,t) ---trans---> (?,?)
        if trans:
            rs = trans(rs, **kwargs)

        return rs
