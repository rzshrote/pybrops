"""
Module implementing selection protocols for weighted genomic selection.
"""

from numbers import Integral
from typing import Callable, Union
import numpy
import types
from pybrops.core.error.error_value_python import check_is_gt

from pybrops.opt.algo.NSGA2SetGeneticAlgorithm import NSGA2SetGeneticAlgorithm
from pybrops.breed.prot.sel.UnconstrainedSelectionProtocol import UnconstrainedSelectionProtocol
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_ndarray
from pybrops.core.error import check_is_callable
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.random.prng import global_prng

class WeightedGenomicSelection(UnconstrainedSelectionProtocol):
    """
    Class implementing selection protocols for weighted genomic selection.

    # TODO: add formulae for methodology.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            objfn_trans = None, 
            objfn_trans_kwargs = None, 
            objfn_wt = 1.0, 
            ndset_trans = None, 
            ndset_trans_kwargs = None, 
            ndset_wt = 1.0, 
            rng = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for WeightedGenomicSelection class.

        Parameters
        ----------
        nparent : int
        ncross : int
        nprogeny : int
        """
        super(WeightedGenomicSelection, self).__init__(**kwargs)

        # variable assignment
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
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
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, method = "single", nparent = None, ncross = None, nprogeny = None, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None, ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = None, **kwargs: dict):
        """
        Select parents individuals for breeding.

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
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        method : str
            Options: "single", "pareto"
        nparent : int, None
            Number of parents. If None, use default.
        ncross : int, None
            Number of crosses per configuration. If None, use default.
        nprogeny : int
            Number of progeny per cross. If None, use default.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing four objects: (pgmat, sel, ncross, nprogeny)

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
        if objfn_trans is None:
            objfn_trans = self.objfn_trans
        if objfn_trans_kwargs is None:
            objfn_trans_kwargs = self.objfn_trans_kwargs
        if objfn_wt is None:
            objfn_wt = self.objfn_wt
        if ndset_trans is None:
            ndset_trans = self.ndset_trans
        if ndset_trans_kwargs is None:
            ndset_trans_kwargs = self.ndset_trans_kwargs
        if ndset_wt is None:
            ndset_wt = self.ndset_wt

        # convert method string to lower
        method = method.lower()

        # single objective method: objfn_trans returns a single value for each
        # selection configuration
        if method == "single":
            # get vectorized objective function
            objfn = self.objfn(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                trans = objfn_trans,
                trans_kwargs = objfn_trans_kwargs
            )

            # get all wGEBVs for each individual
            # (n,)
            wgebv = [objfn(i) for i in range(gmat.ntaxa)]

            # convert to numpy.ndarray
            wgebv = numpy.array(wgebv)

            # multiply the objectives by objfn_wt to transform to maximizing function
            # (n,) * scalar -> (n,)
            wgebv = wgebv * objfn_wt

            # get indices of top nparent GEBVs
            sel = wgebv.argsort()[::-1][:nparent]

            # shuffle indices for random mating
            self.rng.shuffle(sel)

            # get GEBVs for reference
            if miscout is not None:
                miscout["wgebv"] = wgebv

            return pgmat, sel, ncross, nprogeny

        # multi-objective method: objfn_trans returns a multiple values for each
        # selection configuration
        elif method == "pareto":
            # get the pareto frontier
            frontier, sel_config = self.pareto(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                miscout = miscout,
                nparent = nparent,
                objfn_trans = objfn_trans,
                objfn_trans_kwargs = objfn_trans_kwargs,
                objfn_wt = objfn_wt
            )

            # get scores for each of the points along the pareto frontier
            score = ndset_wt * ndset_trans(frontier, **ndset_trans_kwargs)

            # get index of maximum score
            ix = score.argmax()

            # add fields to miscout
            if miscout is not None:
                miscout["frontier"] = frontier
                miscout["sel_config"] = sel_config

            return pgmat, sel_config[ix], ncross, nprogeny
        else:
            raise ValueError("argument 'method' must be either 'single' or 'pareto'")

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, **kwargs: dict):
        """
        Return an objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Not used by this function.
        gmat : GenotypeMatrix
            Used by this function. Input genotype matrix.
        ptdf : PhenotypeDataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Not used by this function.
        gpmod : AdditiveLinearGenomicModel
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

        # get pointers to raw numpy.ndarray matrices
        mat = gmat.mat  # (n,p) get genotype matrix
        u = gpmod.u_a   # (p,t) get regression coefficients

        # calculate weight adjustments for WGS
        afreq = gmat.afreq()[:,None]        # (p,1) allele frequencies
        fafreq = numpy.where(               # (p,t) calculate favorable allele frequencies
            u > 0.0,                        # if dominant (1) allele is beneficial
            afreq,                          # get dominant allele frequency
            1.0 - afreq                     # else get recessive allele frequency
        )
        fafreq[fafreq <= 0.0] = 1.0         # avoid division by zero/imaginary
        uwt = numpy.power(fafreq, -0.5)  # calculate weights: 1/sqrt(p)

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,                 # byte code pointer
            self.objfn_static.__globals__,              # global variables
            None,                                       # new name for the function
            (mat, u, uwt, trans, trans_kwargs),   # default values for arguments
            self.objfn_static.__closure__               # closure byte code pointer
        )

        return outfn

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, **kwargs: dict):
        """
        Return a vectorized objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Not used by this function.
        gmat : GenotypeMatrix
            Used by this function. Input genotype matrix.
        ptdf : PhenotypeDataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Not used by this function.
        gpmod : AdditiveLinearGenomicModel
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

        # get pointers to raw numpy.ndarray matrices
        mat = gmat.mat  # (n,p) get genotype matrix
        u = gpmod.u_a   # (p,t) get regression coefficients

        # calculate weight adjustments for WGS
        afreq = gmat.afreq()[:,None]        # (p,1) allele frequencies
        fafreq = numpy.where(               # (p,t) calculate favorable allele frequencies
            u > 0.0,                        # if dominant (1) allele is beneficial
            afreq,                          # get dominant allele frequency
            1.0 - afreq                     # else get recessive allele frequency
        )
        fafreq[fafreq <= 0.0] = 1.0         # avoid division by zero/imaginary
        uwt = numpy.power(fafreq, -0.5)  # calculate weights: 1/sqrt(p)

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,             # byte code pointer
            self.objfn_vec_static.__globals__,          # global variables
            None,                                       # new name for the function
            (mat, u, uwt, trans, trans_kwargs),   # default values for arguments
            self.objfn_vec_static.__closure__           # closure byte code pointer
        )

        return outfn

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, nparent = None, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None, **kwargs: dict):
        """
        Calculate a Pareto frontier for objectives.

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

        Returns
        -------
        out : tuple
            A tuple containing two objects ``(frontier, sel_config)``.

            Where:

            - frontier is a ``numpy.ndarray`` of shape ``(q,v)`` containing
              Pareto frontier points.
            - sel_config is a ``numpy.ndarray`` of shape ``(q,k)`` containing
              parent selection decisions for each corresponding point in the
              Pareto frontier.

            Where:

            - ``q`` is the number of points in the frontier.
            - ``v`` is the number of objectives for the frontier.
            - ``k`` is the number of search space decision variables.
        """
        if nparent is None:
            nparent = self.nparent
        if objfn_trans is None:
            objfn_trans = self.objfn_trans
        if objfn_trans_kwargs is None:
            objfn_trans_kwargs = self.objfn_trans_kwargs
        if objfn_wt is None:
            objfn_wt = self.objfn_wt

        # get number of taxa
        ntaxa = gmat.ntaxa

        # create objective function
        objfn = self.objfn(
            pgmat = pgmat,
            gmat = gmat,
            ptdf = ptdf,
            bvmat = bvmat,
            gpmod = gpmod,
            t_cur = t_cur,
            t_max = t_max,
            trans = objfn_trans,
            trans_kwargs = objfn_trans_kwargs
        )

        # create optimization algorithm
        moalgo = NSGA2SetGeneticAlgorithm(
            rng = self.rng,
            **kwargs
        )

        # TODO: fixme with miscout dictionary
        frontier, sel_config, misc = moalgo.optimize(
            objfn = objfn,                  # objective function
            k = nparent,                    # vector length to optimize (sspace^k)
            sspace = numpy.arange(ntaxa),   # search space options
            objfn_wt = objfn_wt             # weights to apply to each objective
        )

        if miscout is not None:
            for k,i in misc.pairs():
                miscout[k] = i

        return frontier, sel_config

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, mat, u, uwt, trans, kwargs):
        """
        Score a population of individuals based on Weighted Genomic Selection
        (WGS). Scoring for WGS is defined as the sum of weighted Genomic
        Estimated Breeding Values (wGEBV) for a population.

        WGS selects the ``q`` individuals with the largest GEBVs.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape ``(k,)``

            Where:

            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            If ``sel`` is None, use all individuals.
        mat : numpy.ndarray
            A int8 binary genotype matrix of shape ``(n,p)``.

            Where:

            - ``n`` is the number of individuals.
            - ``p`` is the number of markers.
        u : numpy.ndarray
            A trait prediction coefficients matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.
        uwt : numpy.ndarray
            Multiplicative marker weights matrix to apply to the trait
            prediction coefficients provided of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Trait prediction coefficients (:math:`\\textbf{u}`) are transformed as follows:

            .. math::
                \\textbf{u}_{new} = \\textbf{u} \\bigdot \\textbf{uwt}

            Where:

            - :math:`\\bigdot` is the Hadamard product
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        wgs : numpy.ndarray
            A GEBV matrix of shape ``(k,t)`` if ``objwt`` is ``None``.
            A GEBV matrix of shape ``(k,)`` if ``objwt`` shape is ``(t,)``

            Where:

            - ``k`` is the number of individuals selected.
            - ``t`` is the number of traits.
        """
        # if sel is None, slice all individuals
        if sel is None:
            sel = slice(None)

        # CGS calculation explanation
        # Step 1: (n,p) -> (k,p)            # select individuals
        # Step 2: (k,p) . (p,t) -> (k,t)    # calculate wGEBVs
        # Step 3: (k,t).sum(0) -> (t,)      # sum across all individuals
        wgs = mat[sel,:].dot(u*uwt).sum(0)

        # apply transformations
        # (t,) ---trans---> (?,)
        if trans:
            wgs = trans(wgs, **kwargs)

        return wgs

    @staticmethod
    def objfn_vec_static(sel, mat, u, uwt, trans, kwargs):
        """
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
        Genomic Estimated Breeding Values (GEBV) for a population.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape ``(j,k)``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            If ``sel`` is None, score each individual separately: ``(n,1)``
        mat : numpy.ndarray
            A genotype matrix of shape ``(n,p)``.

            Where:

            - ``n`` is the number of individuals.
            - ``p`` is the number of markers.
        u : numpy.ndarray
            A trait prediction coefficients matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.
        uwt : numpy.ndarray
            Multiplicative marker weights matrix to apply to the trait
            prediction coefficients provided of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Trait prediction coefficients (:math:`\\textbf{u}`) are transformed as follows:

            .. math::
                \\textbf{u}_{new} = \\textbf{u} \\bigdot \\textbf{uwt}

            Where:

            - :math:`\\bigdot` is the Hadamard product
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
            A GEBV matrix of shape ``(j,t)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``t`` is the number of traits.
        """
        # if sel is None, slice all individuals
        if sel is None:
            n = mat.shape[0]
            sel = numpy.arange(n).reshape(n,1)

        # CGS calculation explanation
        # (n,p)[(j,k),:] -> (j,k,p)     # select configurations
        # (j,k,p) . (p,t) -> (j,k,t)    # calculate wGEBVs
        # (j,k,t).sum(1) -> (j,t)       # sum across all individuals in config
        cgs = mat[sel,:].dot(u*uwt).sum(1)

        # apply transformations
        # (j,t) ---trans---> (?,?)
        if trans:
            cgs = trans(cgs, **kwargs)

        return cgs
