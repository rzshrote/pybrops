"""
Module implementing selection protocols for binary weighted general 2-norm genomic selection
"""

import types
import numpy
from typing import Union
from typing import Callable

from pybrops.algo.opt.NSGA2SetGeneticAlgorithm import NSGA2SetGeneticAlgorithm
from pybrops.algo.opt.SteepestAscentSetHillClimber import SteepestAscentSetHillClimber
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.breed.prot.sel.targetfn import target_positive
from pybrops.breed.prot.sel.weightfn import weight_absolute
from pybrops.core.error import check_is_callable
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_str
from pybrops.core.error import check_isinstance
from pybrops.core.error import check_is_gt
from pybrops.core.random.prng import global_prng
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class BinaryWeightedGeneralized2NormGenomicSelection(SelectionProtocol):
    """
    docstring for BinaryWeightedGeneralized2NormGenomicSelection.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            nparent: int,
            ncross: int,
            nprogeny: int,
            target: Union[numpy.ndarray,Callable] = target_positive,
            method: str = "single",
            objfn_trans = None,
            objfn_trans_kwargs = None, 
            objfn_wt = -1.0,
            ndset_trans = None, 
            ndset_trans_kwargs = None, 
            ndset_wt = -1.0,
            rng = global_prng, 
            soalgo = None, 
            moalgo = None,
            **kwargs
        ):
        """
        Constructor for BinaryWeightedGeneralized2NormGenomicSelection.
        
        Parameters
        ----------
        nparent : int
            Number of parents to select.
        """
        super(BinaryWeightedGeneralized2NormGenomicSelection, self).__init__(**kwargs)

        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.target = target
        self.method = method
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = objfn_trans_kwargs
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs
        self.ndset_wt = ndset_wt
        self.rng = rng
        self.soalgo = soalgo
        self.moalgo = moalgo
    
    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def nparent():
        doc = "The nparent property."
        def fget(self):
            return self._nparent
        def fset(self, value):
            check_is_int(value, "nparent")      # must be int
            check_is_gt(value, "nparent", 0)    # int must be >0
            self._nparent = value
        def fdel(self):
            del self._nparent
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    nparent = property(**nparent())

    def ncross():
        doc = "The ncross property."
        def fget(self):
            return self._ncross
        def fset(self, value):
            check_is_int(value, "ncross")       # must be int
            check_is_gt(value, "ncross", 0)     # int must be >0
            self._ncross = value
        def fdel(self):
            del self._ncross
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    ncross = property(**ncross())

    def nprogeny():
        doc = "The nprogeny property."
        def fget(self):
            return self._nprogeny
        def fset(self, value):
            check_is_int(value, "nprogeny")     # must be int
            check_is_gt(value, "nprogeny", 0)   # int must be >0
            self._nprogeny = value
        def fdel(self):
            del self._nprogeny
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    nprogeny = property(**nprogeny())

    def target():
        doc = "The target property."
        def fget(self):
            return self._target
        def fset(self, value):
            check_isinstance(value, "target", (numpy.ndarray, Callable))
            self._target = value
        def fdel(self):
            del self._target
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    target = property(**target())

    def method():
        doc = "The method property."
        def fget(self):
            return self._method
        def fset(self, value):
            check_is_str(value, "method")       # must be string
            value = value.lower()               # convert to lowercase
            options = ("single", "pareto")      # method options
            if value not in options:            # if not method supported
                raise ValueError(               # raise ValueError
                    "Unsupported 'method'. Options are: " +
                    ", ".join(map(str, options))
                )
            self._method = value
        def fdel(self):
            del self._method
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    method = property(**method())

    def objfn_trans():
        doc = "The objfn_trans property."
        def fget(self):
            return self._objfn_trans
        def fset(self, value):
            if value is not None:                       # if given object
                check_is_callable(value, "objfn_trans") # must be callable
            self._objfn_trans = value
        def fdel(self):
            del self._objfn_trans
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    objfn_trans = property(**objfn_trans())

    def objfn_trans_kwargs():
        doc = "The objfn_trans_kwargs property."
        def fget(self):
            return self._objfn_trans_kwargs
        def fset(self, value):
            if value is None:                           # if given None
                value = {}                              # set default to empty dict
            check_is_dict(value, "objfn_trans_kwargs")  # check is dict
            self._objfn_trans_kwargs = value
        def fdel(self):
            del self._objfn_trans_kwargs
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    objfn_trans_kwargs = property(**objfn_trans_kwargs())

    def objfn_wt():
        doc = "The objfn_wt property."
        def fget(self):
            return self._objfn_wt
        def fset(self, value):
            self._objfn_wt = value
        def fdel(self):
            del self._objfn_wt
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    objfn_wt = property(**objfn_wt())

    def ndset_trans():
        doc = "The ndset_trans property."
        def fget(self):
            return self._ndset_trans
        def fset(self, value):
            if value is not None:                       # if given object
                check_is_callable(value, "ndset_trans") # must be callable
            self._ndset_trans = value
        def fdel(self):
            del self._ndset_trans
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    ndset_trans = property(**ndset_trans())

    def ndset_trans_kwargs():
        doc = "The ndset_trans_kwargs property."
        def fget(self):
            return self._ndset_trans_kwargs
        def fset(self, value):
            if value is None:                           # if given None
                value = {}                              # set default to empty dict
            check_is_dict(value, "ndset_trans_kwargs")  # check is dict
            self._ndset_trans_kwargs = value
        def fdel(self):
            del self._ndset_trans_kwargs
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    ndset_trans_kwargs = property(**ndset_trans_kwargs())

    def ndset_wt():
        doc = "The ndset_wt property."
        def fget(self):
            return self._ndset_wt
        def fset(self, value):
            self._ndset_wt = value
        def fdel(self):
            del self._ndset_wt
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    ndset_wt = property(**ndset_wt())

    def rng():
        doc = "The rng property."
        def fget(self):
            return self._rng
        def fset(self, value):
            if value is None:
                value = global_prng
            check_is_Generator_or_RandomState(value, "rng") # check is numpy.Generator
            self._rng = value
        def fdel(self):
            del self._rng
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    rng = property(**rng())

    def soalgo():
        doc = "The soalgo property."
        def fget(self):
            return self._soalgo
        def fset(self, value):
            if value is None:
                value = SteepestAscentSetHillClimber(
                    rng = self.rng  # PRNG source
                )
            self._soalgo = value
        def fdel(self):
            del self._soalgo
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    soalgo = property(**soalgo())

    def moalgo():
        doc = "The moalgo property."
        def fget(self):
            return self._moalgo
        def fset(self, value):
            if value is None:
                value = NSGA2SetGeneticAlgorithm(
                    ngen = 250,     # number of generations to evolve
                    mu = 100,       # number of parents in population
                    lamb = 100,     # number of progeny to produce
                    M = 1.5,        # algorithm crossover genetic map length
                    rng = self.rng  # PRNG source
                )
            self._moalgo = value
        def fdel(self):
            del self._moalgo
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    moalgo = property(**moalgo())

    ############################################################################
    ########################## Private Object Methods ##########################
    ############################################################################
    def _calc_C(self, gmat: GenotypeMatrix, gpmod: AdditiveLinearGenomicModel):
        # declare target allele frequency and marker weight variables
        u_a = gpmod.u_a # (p,t)
        afreq = gmat.afreq()[:,None] # (p,1)
        tfreq = self.target(u_a) if callable(self.target) else self.target # (p,t)

        # calculate marker weights based on:
        # w = |u_a| / sqrt( (p^t) * (1-p)^(1-t) )
        # w = |u_a| * p^(-t/2) * (1-p)^((t-1)/2)
        # u_a = marker effects; p = allele frequencies; t = target allele frequencies

        # calculate the denominator
        denom = numpy.power(afreq, tfreq) * numpy.power(1.0 - afreq, 1.0 - tfreq)
        denom[denom == 0.0] = 1.0 # prevents division by zero
        denom = numpy.power(denom, -0.5) # take square root

        # calculate marker weights
        mkrwt = numpy.absolute(u_a) * denom

        # make sure mkrwt and tfreq have the same dimensions
        if mkrwt.shape != tfreq.shape:
            raise ValueError("marker weights and target allele frequencies do not have the same shape")

        # get number of traits and taxa
        ntrait = mkrwt.shape[1]
        ntaxa = gmat.ntaxa
        ploidy = float(gmat.ploidy)

        # allocate a tensor for storing distance matrices
        C = numpy.empty((ntrait,ntaxa,ntaxa), dtype = "float64")

        # get the genotype matrix as {0,1,2,...}
        # (n,p)
        X = gmat.tacount()

        # calculate a distance matrix for each trait
        for trait in range(ntrait):
            # multiply afreq by ploidy level to get the number of alleles on which to center
            # (p,t)[None,:,ix] -> (1,p)
            # scalar * (1,p) -> (1,p)
            M = ploidy * tfreq[None,:,trait]

            # calculate the Z matrix
            # (n,p) - (1,p) -> (n,p)
            Z = X - M

            # get marker weights
            # (p,t)[None,:,ix] -> (1,p)
            W = mkrwt[None,:,trait]

            # calculate the weighted G matrix for the current trait
            # (n,p) * (1,p) -> (n,p)
            # (n,p) @ (p,n) -> (n,n)
            G = (Z * W).dot(Z.T)

            # apply jitter to G matrix if needed
            diagix = numpy.diag_indices_from(G)
            mat_diag_old = G[diagix].copy()
            counter = 0
            is_not_posdef = not numpy.all(numpy.linalg.eigvals(G) >= 2e-14)

            # attempt to apply jitter in nattempt or less attempts
            while (is_not_posdef) and (counter < 10):
                G[diagix] = mat_diag_old + numpy.random.uniform(1e-10, 1e-6, len(mat_diag_old))
                is_not_posdef = not numpy.all(numpy.linalg.eigvals(G) >= 2e-14)
                counter += 1
            
            # if we still weren't able to find appropriate jitter, then give old diagonals and warn
            if is_not_posdef:
                G[diagix] = mat_diag_old
                raise ValueError("Unable to successfully apply jitter to genomic relationship matrix meet eigenvalue tolerance")

            # take Cholesky decomposition to get upper triangle 
            C[trait,:,:] = numpy.linalg.cholesky(G).T

        # return Cholesky decomposition tensor
        return C

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
        # Solve problem using a single objective method
        if self.method == "single":
            # get vectorized objective function
            objfn = self.objfn(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max
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

            return pgmat, sel, self.ncross, self.nprogeny

        # estimate Pareto frontier, then choose from non-dominated points.
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
                miscout = miscout,
                **kwargs
            )

            # get scores for each of the points along the pareto frontier
            score = self.ndset_wt * self.ndset_trans(frontier, **self.ndset_trans_kwargs)

            # get index of maximum score
            ix = score.argmax()

            # add fields to miscout
            if miscout is not None:
                miscout["frontier"] = frontier
                miscout["sel_config"] = sel_config

            return pgmat, sel_config[ix], self.ncross, self.nprogeny

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return an objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Input phased genotype matrix containing genomes.
        gmat : GenotypeMatrix
            Input genotype matrix.
        ptdf : PhenotypeDataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Input breeding value matrix.
        gpmod : LinearGenomicModel
            Linear genomic prediction model.

        Returns
        -------
        outfn : function
            A selection objective function for the specified problem.
        """
        # get default parameters
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        C = self._calc_C(gmat, gpmod)       # get Cholesky decomposition of genomic relationship matrix: (t,n,n)

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,     # byte code pointer
            self.objfn_static.__globals__,  # global variables
            None,                           # new name for the function
            (C, trans, trans_kwargs),       # default values for last 3 arguments
            self.objfn_static.__closure__   # closure byte code pointer
        )

        return outfn

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return a vectorized objective function.
        """
        # get default parameters
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        C = self._calc_C(gmat, gpmod)       # get Cholesky decomposition of genomic relationship matrix: (t,n,n)

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,     # byte code pointer
            self.objfn_vec_static.__globals__,  # global variables
            None,                               # new name for the function
            (C, trans, trans_kwargs),           # default values for last 3 arguments
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
            t_max = t_max
        )

        # use multi-objective optimization to approximate Pareto front.
        frontier, sel_config, misc = self.moalgo.optimize(
            objfn = objfn,                  # objective function
            k = nparent,                    # vector length to optimize (sspace^k)
            sspace = numpy.arange(ntaxa),   # search space options
            objfn_wt = objfn_wt,            # weights to apply to each objective
            **kwargs
        )

        # handle miscellaneous output
        if miscout is not None:     # if miscout is provided
            miscout.update(misc)    # add 'misc' to 'miscout', overwriting as needed

        return frontier, sel_config

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel: numpy.ndarray, C: numpy.ndarray, trans: Callable, kwargs: dict):
        """
        Score a parent contribution vector according to its distance from a utopian point.
        The goal is to minimize this function.

        Parameters
        ----------
        sel : numpy.ndarray
            A selection indices matrix of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            If ``sel`` is ``None``, use all individuals.
        C : numpy.ndarray
            An upper triangle matrix of shape (n,n) resulting from a Cholesky 
            decomposition of a kinship matrix: K = C'C.

            Where:

            - ``n`` is the number of individuals.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        dist : numpy.ndarray
            A matrix of shape (t,) if ``trans`` is ``None``.

            Distances for each target.

            Where:

            - ``t`` is the number of traits.
        """
        # calculate vector
        # (t,n,n)[:,:,(k,)] -> (t,n,k)
        # (1/k) * (t,n,k).sum(2) -> (t,n,)
        Cx = (1.0 / sel.shape[0]) * C[:,:,sel].sum(2)

        # norm2( (t,n) ) -> (t,)
        dist = numpy.linalg.norm(Cx, ord = 2, axis = 1)

        # apply transformations if needed
        if trans:
            dist = trans(dist, **kwargs)
        
        return dist

    @staticmethod
    def objfn_vec_static(sel: numpy.ndarray, C: numpy.ndarray, trans: Callable, kwargs: dict):
        """
        Score a parent contribution vector according to its distance from a utopian point.
        The goal is to minimize this function.

        Parameters
        ----------
        sel : numpy.ndarray
            A selection indices matrix of shape ``(j,k)``.

            Where:

            - ``j`` is the number of configurations to score.
            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            ``sel`` cannot be ``None``.
        C : numpy.ndarray
            An upper triangle matrix of shape (n,n) resulting from a Cholesky 
            decomposition of a kinship matrix: K = C'C.

            Where:

            - ``n`` is the number of individuals.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        dist : numpy.ndarray
            A matrix of shape (j,t) if ``trans`` is ``None``.

            Distances for each target.

            Where:

            - ``t`` is the number of traits.
        """
        # calculate MEH
        # (t,n,n)[:,:,(j,k)] -> (t,n,j,k)
        # (t,n,j,k).sum(2) -> (t,n,j)
        Cx = (1.0 / sel.shape[1]) * C[:,:,sel].sum(3)
        
        # norm2( (t,n,j), axis=1 ) -> (t,j)
        # (t,j).T -> (t,j)
        dist = numpy.linalg.norm(Cx, ord = 2, axis = 0).T

        # apply transformations if needed
        if trans:
            dist = trans(dist, **kwargs)
        
        return dist