"""
Module implementing selection protocols for optimal haploid value selection.
"""

import types
from typing import Callable, Union
import numpy

from pybrops.algo.opt.NSGA2SetGeneticAlgorithm import NSGA2SetGeneticAlgorithm
from pybrops.algo.opt.OptimizationAlgorithm import OptimizationAlgorithm, check_is_OptimizationAlgorithm
from pybrops.algo.opt.SteepestAscentSetHillClimber import SteepestAscentSetHillClimber
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.core.error import check_is_bool
from pybrops.core.error import check_is_callable
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_gt
from pybrops.core.error import check_is_str
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.random.prng import global_prng
from pybrops.core.util.arrayix import triuix
from pybrops.core.util.arrayix import triudix
from pybrops.core.util.haplo import calc_haplobin
from pybrops.core.util.haplo import calc_haplobin_bounds
from pybrops.core.util.haplo import calc_nhaploblk_chrom

class OptimalHaploidValueSelection(SelectionProtocol):
    """
    Class implementing selection protocols for optimal haploid value selection.

    # TODO: add formulae for methodology.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            nconfig: int, 
            nparent: int, 
            ncross: int, 
            nprogeny: int, 
            nhaploblk: int,
            unique_parents = True, 
            method = "single",
            objfn_trans = None, 
            objfn_trans_kwargs = None, 
            objfn_wt = 1.0,
            ndset_trans = None, 
            ndset_trans_kwargs = None, 
            ndset_wt = 1.0,
            rng = None, 
            soalgo = None, 
            moalgo = None,
            **kwargs : dict
        ):
        """
        Constructor for Optimal Haploid Value Selection (OHV).

        Parameters
        ----------
        nconfig : int
            Number of cross configurations to consider.

            Examples:

            - 20 two-way crosses would be: ``nconfig = 20``
            - 20 three way crosses would be: ``nconfig = 20``
        nparent : int
            Number of parents to per configuration.

            Example:

            - 20 two-way crosses would be: ``nparent = 2``
            - 20 three-way crosses would be: ``nparent = 3``
        ncross : int
            Number of crosses to perform per configuration.
        nprogeny : int
            Number of progeny to derive from each cross configuration.
        nhaploblk : int
            Number of haplotype blocks to segment the genome into.
        unique_parents : bool, default = True
            Whether to allow force unique parents or not.
            If ``True``, all parents in the mating configuration must be unique.
            If ``False``, non-unique parents are allowed. In this scenario,
            self-fertilization is considered as a viable option.
        method : str
            Method of selecting parents.

            +--------------+---------------------------------------------------+
            | Method       | Description                                       |
            +==============+===================================================+
            | ``"single"`` | OHV is transformed to a single objective and      |
            |              | optimization is done on the transformed function. |
            |              | This is done using the ``trans`` function         |
            |              | provided::                                        |
            |              |                                                   |
            |              |    optimize : objfn_trans(MOGS)                   |
            +--------------+---------------------------------------------------+
            | ``"pareto"`` | OHV is transformed by a transformation function,  |
            |              | but NOT reduced to a single objective. The Pareto |
            |              | frontier for this transformed function is mapped  |
            |              | using a multi-objective genetic algorithm.        |
            |              |                                                   |
            |              | Objectives are scaled to :math:`[0,1]` and a      |
            |              | vector orthogonal to the hyperplane defined by    |
            |              | the extremes of the front is drawn starting at    |
            |              | the point defined by ``ndset_trans``. The closest |
            |              | point on the Pareto frontier to the orthogonal    |
            |              | vector is selected.                               |
            +--------------+---------------------------------------------------+
        objfn_trans : function, callable, None
        objfn_trans_kwargs : dict, None
        objfn_wt : float, numpy.ndarray
        ndset_trans : function, callable, None
        ndset_trans_kwargs : dict, None
        ndset_wt : float
        soalgo : OptimizationAlgorithm
            Single-objective optimization algorithm to optimize the objective
            function. If ``None``, use a SteepestAscentSetHillClimber with the
            following parameters::

                soalgo = SteepestAscentSetHillClimber(
                    rng = self.rng  # PRNG source
                )
        moalgo : OptimizationAlgorithm
            Multi-objective optimization algorithm to optimize the objective
            functions. If ``None``, use a NSGA2SetGeneticAlgorithm with the
            following parameters::

                moalgo = NSGA2SetGeneticAlgorithm(
                    ngen = 250,     # number of generations to evolve
                    mu = 100,       # number of parents in population
                    lamb = 100,     # number of progeny to produce
                    M = 1.5,        # algorithm crossover genetic map length
                    rng = self.rng  # PRNG source
                )
        rng : numpy.random.Generator, numpy.random.RandomState
        """
        super(OptimalHaploidValueSelection, self).__init__(**kwargs)

        # error checks and assignments (ORDER DEPENDENT!!!)
        self.nconfig = nconfig
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.nhaploblk = nhaploblk
        self.unique_parents = unique_parents
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
    def nconfig(self) -> int:
        """Number of cross configurations to consider."""
        return self._nconfig
    @nconfig.setter
    def nconfig(self, value: int) -> None:
        """Set number of cross configurations to consider."""
        check_is_int(value, "nconfig")      # must be int
        check_is_gt(value, "nconfig", 0)    # int must be >0
        self._nconfig = value
    @nconfig.deleter
    def nconfig(self) -> None:
        """Delete number of cross configurations to consider."""
        del self._nconfig

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
    def unique_parents(self) -> bool:
        """Whether parents should be unique."""
        return self._unique_parents
    @unique_parents.setter
    def unique_parents(self, value: bool) -> None:
        """Set whether parents should be unique."""
        check_is_bool(value, "unique_parents")
        self._unique_parents = value
    @unique_parents.deleter
    def unique_parents(self) -> None:
        """Delete whether parents should be unique."""
        del self._unique_parents

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
        nblk = calc_nhaploblk_chrom(self.nhaploblk, genpos, chrgrp_stix, chrgrp_spix)

        # ensure there are enough markers per chromosome
        if numpy.any(nblk > chrgrp_len):
            raise RuntimeError(
                "number of haplotype blocks assigned to a chromosome greater than number of available markers"
            )

        # calculate haplotype bins
        hbin = calc_haplobin(nblk, genpos, chrgrp_stix, chrgrp_spix)

        # define shape
        # (m,n,b,t)
        s = (mat.shape[0], mat.shape[1], self.nhaploblk, u.shape[1])

        # allocate haplotype matrix
        # (m,n,b,t)
        hmat = numpy.empty(s, dtype = u.dtype)

        # get boundary indices
        hstix, hspix, hlen = calc_haplobin_bounds(hbin)

        # OPTIMIZE: perhaps eliminate one loop using dot function
        # fill haplotype matrix
        for i in range(hmat.shape[3]):                          # for each trait
            for j,(st,sp) in enumerate(zip(hstix,hspix)):       # for each haplotype block
                hmat[:,:,j,i] = mat[:,:,st:sp].dot(u[st:sp,i])  # take dot product and fill

        return hmat

    def _calc_xmap(self, ntaxa):
        """
        Calculate the cross map.

        Parameters
        ----------
        ntaxa : int
            Number of taxa.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(s,d)`` containing cross map indices.

            Where:

            - ``s`` is the number of elements in the upper triangle, including
              or not including the diagonal (depending on ``unique_parents``).
            - ``d`` is the number of parents in the cross.
        """
        if self.unique_parents:         # if we want unique parents
            return numpy.array(         # create a numpy.ndarray
                list(                   # convert to list
                    triudix(            # generator for indices without diagonal
                        ntaxa,          # number of taxa
                        self.nparent    # number of parents
                    )
                )
            )
        else:                           # otherwise we don't want unique parents
            return numpy.array(         # create a numpy.ndarray
                list(                   # convert to list
                    triuix(             # generator for indices with diagonal
                        ntaxa,          # number of taxa
                        self.nparent    # number of parents
                    )
                )
            )

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs: dict):
        """
        Select parents individuals for breeding.

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
        # single-objective method: objfn_trans returns a single value for each
        # selection configuration
        if self.method == "single":
            # get vectorized objective function
            objfn = self.objfn(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                **kwargs
            )

            # calculate xmap
            xmap = self._calc_xmap(pgmat.ntaxa)

            # get all OHVs for each configuration
            # (s,)
            ohv = numpy.array([objfn([i]) for i in range(len(xmap))])

            # multiply the objectives by objfn_wt to transform to maximizing function
            # (n,) * scalar -> (n,)
            ohv = ohv * self.objfn_wt

            # get indices of top nconfig OHVs
            sel = ohv.argsort()[::-1][:self.nconfig]

            # shuffle indices for random mating
            self.rng.shuffle(sel)

            # convert 'sel' to parent selections (ordered)
            # (kd,)
            sel = xmap[sel,:].flatten()

            # add optimization details to miscellaneous output if miscout was provided
            if miscout is not None:
                miscout["ohv"] = ohv

            return pgmat, sel, self.ncross, self.nprogeny

        # multi-objective method: objfn_trans returns a multiple values for each
        # selection configuration
        elif self.method == "pareto":
            # get the pareto frontier
            frontier, sel_config = self.pareto(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                **kwargs
            )

            # get scores for each of the points along the pareto frontier
            score = self.ndset_wt * self.ndset_trans(frontier, **self.ndset_trans_kwargs)

            # get index of maximum score
            ix = score.argmax()

            # get selection
            sel = sel_config[ix]

            # calculate xmap
            xmap = self._calc_xmap(pgmat.ntaxa)

            # convert 'sel' to parent selections (ordered)
            # (kd,)
            sel = xmap[sel,:].flatten()

            # add optimization details to miscellaneous output if miscout was provided
            if miscout is not None:
                miscout["frontier"] = frontier
                miscout["sel_config"] = sel_config

            return pgmat, sel, self.ncross, self.nprogeny

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return a parent selection objective function.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Phased genotype matrix.
        gmat : GenotypeMatrix
            Not used by this function.
        ptdf : PhenotypeDataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Not used by this function.
        gpmod : LinearGenomicModel
            Linear genomic prediction model.
        """
        ntaxa = pgmat.ntaxa                     # get number of taxa
        xmap = self._calc_xmap(ntaxa)           # (s,p) get the cross map
        mat = self._calc_hmat(pgmat, gpmod)     # (m,n,b,t) get haplotype matrix
        ploidy = pgmat.ploidy                   # (scalar) get ploidy
        trans = self.objfn_trans                # get transformation function
        trans_kwargs = self.objfn_trans_kwargs  # get transformation function keyword arguments

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,                 # byte code pointer
            self.objfn_static.__globals__,              # global variables
            None,                                       # new name for the function
            (xmap, mat, ploidy, trans, trans_kwargs),   # default values for arguments
            self.objfn_static.__closure__               # closure byte code pointer
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
        xmap = self._calc_xmap()                # (s,p) get the cross map
        ploidy = pgmat.ploidy                   # (scalar) get ploidy
        trans = self.objfn_trans                # get transformation function
        trans_kwargs = self.objfn_trans_kwargs  # get transformation function keyword arguments

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,             # byte code pointer
            self.objfn_vec_static.__globals__,          # global variables
            None,                                       # new name for the function
            (xmap, mat, ploidy, trans, trans_kwargs),   # default values for arguments
            self.objfn_vec_static.__closure__           # closure byte code pointer
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
        # create objective function
        objfn = self.objfn(
            pgmat = pgmat,
            gmat = gmat,
            ptdf = ptdf,
            bvmat = bvmat,
            gpmod = gpmod,
            t_cur = t_cur,
            t_max = t_max,
            **kwargs
        )

        # use multi-objective optimization to approximate Pareto front.
        frontier, sel_config, misc = self.moalgo.optimize(
            objfn = objfn,                      # objective function
            k = self.nparent,                   # vector length to optimize (sspace^k)
            sspace = numpy.arange(pgmat.ntaxa), # search space options
            objfn_wt = self.objfn_wt            # weights to apply to each objective
        )

        # handle miscellaneous output
        if miscout is not None:     # if miscout is provided
            miscout.update(misc)    # add 'misc' to 'miscout', overwriting as needed

        # TODO: fix sel_config output format: currently in xmap indices format.
        return frontier, sel_config

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, xmap, mat, ploidy, trans, kwargs):
        """
        Score a population of individuals based on Optimal Haploid Value
        Selection (OHV). Scoring for OHV is defined as the sum of maximum
        haploid breeding values obtainable from a population.

        OHV selects the ``q`` individuals with the largest OHVs.

        Parameters
        ----------
        sel : numpy.ndarray
            A cross selection indices array of shape ``(k,)``.

            Where:

            - ``k`` is the number of crosses to select.
        xmap : numpy.ndarray
            A cross selection index map array of shape ``(s,p)``.

            Where:

            - ``s`` is the size of the sample space (number of cross
              combinations for ``d`` parents).
            - ``p`` is the number of parents.
        mat : numpy.ndarray
            A haplotype effect matrix of shape ``(m,n,b,t)``.

            Where:

            - ``m`` is the number of chromosome phases (2 for diploid, etc.).
            - ``n`` is the number of individuals.
            - ``b`` is the number of haplotype blocks.
            - ``t`` is the number of traits.
        ploidy : int
            Ploidy level of the species.
            In many cases, this should be equal to ``m`` from the ``mat``
            parameter. In cases where data is unphased (``m == 1``), then this
            parameter should be different from ``m``.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard::

                trans(numpy.ndarray, **kwargs: dict):
                    return (scalar or numpy.ndarray)
        kwargs : dict
            Dictionary of keyword arguments to pass to the ``trans`` function.

        Returns
        -------
        ohv : numpy.ndarray
            A OHV matrix of shape ``(t,)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.
        """
        # get the cross configurations
        # (s,d)[(k,),:] -> (k,d)
        sel = xmap[sel,:]

        # get maximum haplotype value
        # multiply by ploidy
        # (m,n,b,t)[:,(k,d),:,:] -> (m,k,d,b,t) # select k individuals
        # (m,k,d,b,t).max((0,2)) -> (k,b,t)     # find maximum haplotype across all parental phases
        # (k,b,t).sum((0,1)) -> (t,)            # add maximum haplotypes for k crosses and b blocks
        # scalar * (t,) -> (t,)                 # multiply result by number of phases
        ohv = ploidy * mat[:,sel,:,:].max((0,2)).sum((0,1))

        # apply transformations
        # (t,) ---trans---> (?,)
        if trans:
            ohv = trans(ohv, **kwargs)

        return ohv

    @staticmethod
    def objfn_vec_static(sel, xmap, mat, ploidy, trans, kwargs):
        """
        Score a population of individuals based on Optimal Haploid Value
        Selection (OHV). Scoring for OHV is defined as the sum of maximum
        haploid breeding values obtainable from a population.

        OHV selects the ``q`` individuals with the largest OHVs.

        Parameters
        ----------
        sel : numpy.ndarray
            A selection indices matrix of shape ``(j,k)``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``k`` is the number of individuals to select.
        xmap : numpy.ndarray
            A cross selection index map array of shape ``(s,p)``.

            Where:

            - ``s`` is the size of the sample space (number of cross
              combinations for ``d`` parents).
            - ``p`` is the number of parents.
        mat : numpy.ndarray
            A haplotype effect matrix of shape ``(t,m,n,b)``.

            Where:

            - ``t`` is the number of traits.
            - ``m`` is the number of chromosome phases (2 for diploid, etc.).
            - ``n`` is the number of individuals.
            - ``b`` is the number of haplotype blocks.
        ploidy : int
            Ploidy level of the species.
            In many cases, this should be equal to ``m`` from the ``mat``
            parameter. In cases where data is unphased (``m == 1``), then this
            parameter should be different from ``m``.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard::

                trans(numpy.ndarray, **kwargs: dict):
                    return (scalar or numpy.ndarray)
        kwargs : dict
            Dictionary of keyword arguments to pass to the ``trans`` function.

        Returns
        -------
        ohv : numpy.ndarray
            A OHV matrix of shape ``(t,)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.
        """
        # get the cross configurations
        # (s,d)[(j,k),:] -> (j,k,d)
        sel = xmap[sel,:]

        # get maximum haplotype value
        # multiply by ploidy
        # (m,n,b,t)[:,(j,k,d),:] -> (m,j,k,d,b,t)   # select k individuals
        # (m,j,k,d,b,t).max((0,3)) -> (j,k,b,t)     # find maximum haplotype across all parental phases
        # (j,k,b,t).sum((1,2)) -> (j,t)             # add maximum haplotypes for k crosses and b blocks
        # scalar * (j,t) -> (j,t)                   # multiply result by number of phases
        ohv = ploidy * mat[:,sel,:,:].max((0,3)).sum((1,2))

        # apply transformations
        # (j,t) ---trans---> (?,?)
        if trans:
            ohv = trans(ohv, **kwargs)

        return ohv
