"""
Module implementing selection protocols for Genotype Building selection.
"""

from typing import Callable, Union
import numpy
import types

from pybrops.opt.algo.NSGA2SetGeneticAlgorithm import NSGA2SetGeneticAlgorithm
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm, check_is_OptimizationAlgorithm
from pybrops.opt.algo.SteepestAscentSetHillClimber import SteepestAscentSetHillClimber
from pybrops.breed.prot.sel.UnconstrainedSelectionProtocol import UnconstrainedSelectionProtocol
from pybrops.core.error import check_is_callable
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_gt
from pybrops.core.error import check_is_str
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.error.error_value_python import check_is_lt
from pybrops.core.random.prng import global_prng
from pybrops.core.util.haplo import nhaploblk_chrom
from pybrops.core.util.haplo import haplobin
from pybrops.core.util.haplo import haplobin_bounds

class GenotypeBuildingSelection(UnconstrainedSelectionProtocol):
    """
    Class implementing selection protocols for Genotype Building selection.

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
            nhaploblk: int,
            nbestfndr: int,
            method = "single",
            objfn_trans = None, 
            objfn_trans_kwargs = None, 
            objfn_wt = 1.0,
            ndset_trans = None, 
            ndset_trans_kwargs = None, 
            ndset_wt = 1.0,
            soalgo = None, 
            moalgo = None,
            rng = None, 
            **kwargs: dict
        ):
        """
        Constructor for Genotype Building selection (OPV).

        Parameters
        ----------
        nparent : int
            Number of parents to select.
        ncross : int
            Number of crosses per configuration.
        nprogeny : int
            Number of progeny to derive from each cross.
        nhaploblk : int
            Number of haplotype blocks to divide the genome into.
        objfn_trans : function, callable, None
        objfn_trans_kwargs : dict, None
        objfn_wt : float, numpy.ndarray
        ndset_trans : function, callable, None
        ndset_trans_kwargs : dict, None
        ndset_wt : float
        rng : numpy.random.Generator, numpy.random.RandomState
        """
        super(GenotypeBuildingSelection, self).__init__(**kwargs)

        # error checks and assignments (ORDER DEPENDENT!!!)
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.nhaploblk = nhaploblk
        self.nbestfndr = nbestfndr # must be less than or equal to nparent
        self.method = method
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = objfn_trans_kwargs # property replaces None with {}
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs # property replaces None with {}
        self.ndset_wt = ndset_wt
        self.rng = rng  # property replaces None with pybrops.core.random
        # soalgo, moalgo MUST GO AFTER 'rng'; properties provide default if None
        self.soalgo = soalgo
        self.moalgo = moalgo

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
    def nhaploblk(self) -> int:
        """Number of haplotype blocks to consider."""
        return self._nhaploblk
    @nhaploblk.setter
    def nhaploblk(self, value: int) -> None:
        """Set number of haplotype blocks to consider."""
        check_is_int(value, "nhaploblk")    # must be int
        check_is_gt(value, "nhaploblk", 0)  # int must be >0
        self._nhaploblk = value
    @nhaploblk.deleter
    def nhaploblk(self) -> None:
        """Delete number of haplotype blocks to consider."""
        del self._nhaploblk

    @property
    def nbestfndr(self) -> int:
        """Number of best founders to consider."""
        return self._nbestfndr
    @nbestfndr.setter
    def nbestfndr(self, value: int) -> None:
        """Set number of best founders to consider."""
        check_is_int(value, "nbestfndr")                # must be int
        check_is_gt(value, "nbestfndr", 0)              # int must be >0
        check_is_lt(value, "nbestfndr", self.nparent)   # cannot have more founders than number of parents selected
        self._nbestfndr = value
    @nbestfndr.deleter
    def nbestfndr(self) -> None:
        """Delete number of best founders to consider."""
        del self._nbestfndr

    @property
    def method(self) -> str:
        """Selection method."""
        return self._method
    @method.setter
    def method(self, value: str) -> None:
        """Set selection method."""
        check_is_str(value, "method")       # must be string
        value = value.lower()               # convert to lowercase
        options = ("single", "pareto")      # method options
        # if not method supported raise ValueError
        if value not in options:
            raise ValueError("Unsupported 'method'. Options are: " + ", ".join(map(str, options)))
        self._method = value
    @method.deleter
    def method(self) -> None:
        """Delete selection method."""
        del self._method

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
    def soalgo(self) -> OptimizationAlgorithm:
        """Single objective optimization algorithm."""
        return self._soalgo
    @soalgo.setter
    def soalgo(self, value: Union[OptimizationAlgorithm,None]) -> None:
        """Set single objective optimization algorithm."""
        if value is None:
            value = SteepestAscentSetHillClimber(rng = self.rng)
        check_is_OptimizationAlgorithm(value, "soalgo")
        self._soalgo = value
    @soalgo.deleter
    def soalgo(self) -> None:
        """Delete single objective optimization algorithm."""
        del self._soalgo

    @property
    def moalgo(self) -> OptimizationAlgorithm:
        """Multi-objective opimization algorithm."""
        return self._moalgo
    @moalgo.setter
    def moalgo(self, value: Union[OptimizationAlgorithm,None]) -> None:
        """Set multi-objective opimization algorithm."""
        if value is None:
            value = NSGA2SetGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                mu = 100,       # number of parents in population
                lamb = 100,     # number of progeny to produce
                M = 1.5,        # algorithm crossover genetic map length
                rng = self.rng  # PRNG source
            )
        check_is_OptimizationAlgorithm(value, "moalgo")
        self._moalgo = value
    @moalgo.deleter
    def moalgo(self) -> None:
        """Delete multi-objective opimization algorithm."""
        del self._moalgo

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
    ########################## Private Object Methods ##########################
    ############################################################################
    def _calc_hmat(self, gmat, mod):
        """
        Calculate a haplotype matrix from a genome matrix and model.

        Parameters
        ----------
        gmat : PhasedGenotypeMatrix
            A genome matrix.
        mod : DenseAdditiveLinearGenomicModel
            A genomic prediction model.

        Returns
        -------
        hmat : numpy.ndarray
            A haplotype effect matrix of shape ``(m,n,b,t)``.
        """
        mat         = gmat.mat              # get genotypes
        genpos      = gmat.vrnt_genpos      # get genetic positions
        chrgrp_stix = gmat.vrnt_chrgrp_stix # get chromosome start indices
        chrgrp_spix = gmat.vrnt_chrgrp_spix # get chromosome stop indices
        chrgrp_len  = gmat.vrnt_chrgrp_len  # get chromosome marker lengths
        u           = mod.u_a               # get regression coefficients

        if (chrgrp_stix is None) or (chrgrp_spix is None):
            raise RuntimeError("markers are not sorted by chromosome position")

        # get number of chromosomes
        nchr = len(chrgrp_stix)

        if self.nhaploblk < nchr:
            raise RuntimeError("number of haplotype blocks is less than the number of chromosomes")

        # calculate number of marker blocks to assign to each chromosome
        nblk = nhaploblk_chrom(self.nhaploblk, genpos, chrgrp_stix, chrgrp_spix)

        # ensure there are enough markers per chromosome
        if numpy.any(nblk > chrgrp_len):
            raise RuntimeError(
                "number of haplotype blocks assigned to a chromosome greater than number of available markers"
            )

        # calculate haplotype bins
        hbin = haplobin(nblk, genpos, chrgrp_stix, chrgrp_spix)

        # define shape
        s = (mat.shape[0], mat.shape[1], self.nhaploblk, u.shape[1]) # (m,n,b,t)

        # allocate haplotype matrix
        hmat = numpy.empty(s, dtype = u.dtype)   # (m,n,b,t)

        # get boundary indices
        hstix, hspix, hlen = haplobin_bounds(hbin)

        # OPTIMIZE: perhaps eliminate one loop using dot function
        # fill haplotype matrix
        for i in range(hmat.shape[3]):                          # for each trait
            for j,(st,sp) in enumerate(zip(hstix,hspix)):       # for each haplotype block
                hmat[:,:,j,i] = mat[:,:,st:sp].dot(u[st:sp,i])  # take dot product and fill

        return hmat

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs: dict):
        """
        Select individuals for breeding.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes (unphased most likely)
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
        nparent : int
        ncross : int
        nprogeny : int
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
        # get selection parameters
        nparent = self.nparent
        ncross = self.ncross
        nprogeny = self.nprogeny
        objfn_trans = self.objfn_trans
        objfn_trans_kwargs = self.objfn_trans_kwargs
        objfn_wt = self.objfn_wt
        ndset_trans = self.ndset_trans
        ndset_trans_kwargs = self.ndset_trans_kwargs
        ndset_wt = self.ndset_wt
        method = self.method

        # single-objective method: objfn_trans returns a single value for each
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

            # optimize using single objective algorithm
            sel_score, sel, misc = self.soalgo.optimize(
                objfn,                              # objective function
                k = self.nparent,                   # number of parents to select
                sspace = numpy.arange(pgmat.ntaxa), # parental indices
                objfn_wt = self.objfn_wt,           # maximizing function
                **kwargs
            )

            # shuffle selection to ensure random mating
            numpy.random.shuffle(sel)

            # add optimization details to miscellaneous output
            if miscout is not None:
                miscout["sel_score"] = sel_score
                miscout["sel"] = sel
                miscout.update(misc) # add dict to dict

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

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return a selection objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Input genome matrix.
        gmat : GenotypeMatrix
            Not used by this function.
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
        # get selection parameters
        mat = self._calc_hmat(pgmat, gpmod)     # (m,n,b,t) get haplotype matrix
        nbestfndr = self.nbestfndr              # get number of best founders
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,         # byte code pointer
            self.objfn_static.__globals__,      # global variables
            None,                               # new name for the function
            (mat, nbestfndr, trans, trans_kwargs), # default values for arguments
            self.objfn_static.__closure__       # closure byte code pointer
        )

        return outfn

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return a vectorized selection objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Input genome matrix.
        gmat : GenotypeMatrix
            Not used by this function.
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
        # get selection parameters
        mat = self._calc_hmat(pgmat, gpmod)     # (m,n,b,t) get haplotype matrix
        nbestfndr = self.nbestfndr              # get number of best founders
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,     # byte code pointer
            self.objfn_vec_static.__globals__,  # global variables
            None,                               # new name for the function
            (mat, nbestfndr, trans, trans_kwargs), # default values for arguments
            self.objfn_vec_static.__closure__   # closure byte code pointer
        )

        return outfn

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs: dict):
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

            - ``frontier`` is a ``numpy.ndarray`` of shape ``(q,v)`` containing
              Pareto frontier points.
            - ``sel_config`` is a ``numpy.ndarray`` of shape ``(q,k)`` containing
              parent selection decisions for each corresponding point in the
              Pareto frontier.

            Where:

            - ``q`` is the number of points in the frontier.
            - ``v`` is the number of objectives for the frontier.
            - ``k`` is the number of search space decision variables.
        """
        # get selection parameters
        nparent = self.nparent
        objfn_trans = self.objfn_trans
        objfn_trans_kwargs = self.objfn_trans_kwargs
        objfn_wt = self.objfn_wt

        # get number of taxa
        ntaxa = pgmat.ntaxa

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

        # use multi-objective optimization to approximate Pareto front.
        frontier, sel_config, misc = self.moalgo.optimize(
            objfn = objfn,                  # objective function
            k = nparent,                    # vector length to optimize (sspace^k)
            sspace = numpy.arange(ntaxa),   # search space options
            objfn_wt = objfn_wt             # weights to apply to each objective
        )

        # handle miscellaneous output
        if miscout is not None:     # if miscout is provided
            miscout.update(misc)    # add 'misc' to 'miscout', overwriting as needed

        return frontier, sel_config

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, mat, nbestfndr, trans, kwargs):
        """
        Score a population of individuals based on Genotype Building
        Selection.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            If ``sel`` is None, use all individuals.
        mat : numpy.ndarray
            A haplotype effect matrix of shape ``(m,n,b,t)``.

            Where:

            - ``m`` is the number of chromosome phases (2 for diploid, etc.).
            - ``n`` is the number of individuals.
            - ``h`` is the number of haplotype blocks.
            - ``t`` is the number of traits.
        nbestfndr : int
            Number of best founders to select. Must be less than or equal to ``n``.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        gb : numpy.ndarray
            An GB matrix of shape ``(t,)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.
        """
        # get number of phases
        nphase = mat.shape[0]

        # get best haplotype within each individual
        # (m,n,h,t)[:,(k,),:,:] -> (m,k,h,t)
        # (m,k,h,t).max(0) -> (k,h,t)
        bestphase = mat[:,sel,:,:].max(0)

        # sort the bestphase tensor along the individual axis to get best individuals
        # (k,h,t) 
        bestphase.sort(0)

        # get top haplotype index boundaries
        k = len(sel)
        st = k - nbestfndr
        sp = k

        # get top haplotypes and sum across the individual and haplotype axes
        # scale by nphase / nbestfndr
        # (k,h,t)[(nbestfndr,),:,:] -> (nbestfndr,h,t)
        # (nbestfndr,h,t).sum((0,1)) -> (t,)
        # scalar * (t,) -> (t,)
        gb = (nphase / nbestfndr) * bestphase[st:sp,:,:].sum((0,1))

        # apply transformations
        # (t,) ---trans---> (?,)
        if trans is not None:
            gb = trans(gb, **kwargs)

        return gb

    @staticmethod
    def objfn_vec_static(sel, mat, nbestfndr, trans, kwargs):
        """
        Score a population of individuals based on Genotype Building
        Selection.

        Parameters
        ----------
        sel : numpy.ndarray
            A selection indices matrix of shape ``(j,k)``.

            Where:

            - ``j`` is the number of configurations to score.
            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
        mat : numpy.ndarray
            A haplotype effect matrix of shape ``(m,n,h,t)``.

            Where:

            - ``m`` is the number of chromosome phases (2 for diploid, etc.).
            - ``n`` is the number of individuals.
            - ``h`` is the number of haplotype blocks.
            - ``t`` is the number of traits.
        nbestfndr : int
            Number of best founders to select. Must be less than or equal to ``n``.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        opv : numpy.ndarray
            An OPV matrix of shape ``(j,t)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``t`` is the number of traits.
        """
        # get number of phases
        nphase = mat.shape[0]

        # get best haplotype within each individual
        # (m,n,h,t)[:,(j,k),:,:] -> (m,j,k,h,t)
        # (m,j,k,h,t).max(0) -> (j,k,h,t)
        bestphase = mat[:,sel,:,:].max(0)

        # sort the bestphase tensor along the individual axis to get best individuals
        # (j,k,h,t) 
        bestphase.sort(1)

        # get top haplotype index boundaries
        k = len(sel)
        st = k - nbestfndr
        sp = k

        # get top haplotypes and sum across the individual and haplotype axes
        # scale by nphase / nbestfndr
        # (j,k,h,t)[:,(nbestfndr,),:,:] -> (j,nbestfndr,h,t)
        # (j,nbestfndr,h,t).sum((1,2)) -> (j,t)
        # scalar * (j,t) -> (j,t)
        gb = (nphase / nbestfndr) * bestphase[:,st:sp,:,:].sum((1,2))

        # apply transformations
        # (j,t) ---trans---> (j,?)
        if trans is not None:
            gb = trans(gb, **kwargs)

        return gb
