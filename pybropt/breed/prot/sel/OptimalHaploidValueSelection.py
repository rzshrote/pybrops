import numpy
from sklearn.cluster import KMeans

from pybropt.breed.prot.sel.SelectionProtocol import SelectionProtocol

import pybropt.core.random
from pybropt.core.error import check_is_int
from pybropt.core.error import check_is_Generator
from pybropt.core.error import check_is_bool
from pybropt.core.error import check_is_gt
from pybropt.core.error import check_is_callable
from pybropt.core.error import check_is_dict
from pybropt.core.util.arrayix import triuix
from pybropt.core.util.arrayix import triudix
from pybropt.core.util.haplo import calc_nhaploblk_chrom
from pybropt.core.util.haplo import calc_haplobin
from pybropt.core.util.haplo import calc_haplobin_bounds

class OptimalHaploidValueSelection(SelectionProtocol):
    """docstring for OptimalHaploidValueSelection."""

    def __init__(self,
    nconfig, nparent, ncross, nprogeny, nhaploblk, unique_parents = True,
    objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = 1.0,
    ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = 1.0,
    rng = None, **kwargs):
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
        objfn_trans : function, callable, None
        objfn_trans_kwargs : dict, None
        objfn_wt : float, numpy.ndarray
        ndset_trans : function, callable, None
        ndset_trans_kwargs : dict, None
        ndset_wt : float
        rng : numpy.Generator
        """
        super(OptimalHaploidValueSelection, self).__init__(**kwargs)

        # error checks and assignments
        self.nconfig = nconfig
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.nhaploblk = nhaploblk
        self.unique_parents = unique_parents
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = objfn_trans_kwargs # property replaces None with {}
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs # property replaces None with {}
        self.ndset_wt = ndset_wt
        self.rng = rng  # property replaces None with pybropt.core.random

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def nconfig():
        doc = "The nconfig property."
        def fget(self):
            return self._nconfig
        def fset(self, value):
            check_is_int(value, "nconfig")      # must be int
            check_is_gt(value, "nconfig", 0)    # int must be >0
            self._nconfig = value
        def fdel(self):
            del self._nconfig
        return locals()
    nconfig = property(**nconfig())

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
        return locals()
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
        return locals()
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
        return locals()
    nprogeny = property(**nprogeny())

    def nhaploblk():
        doc = "The nhaploblk property."
        def fget(self):
            return self._nhaploblk
        def fset(self, value):
            check_is_int(value, "nhaploblk")    # must be int
            check_is_gt(value, "nhaploblk", 0)  # int must be >0
            self._nhaploblk = value
        def fdel(self):
            del self._nhaploblk
        return locals()
    nhaploblk = property(**nhaploblk())

    def unique_parents():
        doc = "The unique_parents property."
        def fget(self):
            return self._unique_parents
        def fset(self, value):
            check_is_bool(value, "unique_parents")
            self._unique_parents = value
        def fdel(self):
            del self._unique_parents
        return locals()
    unique_parents = property(**unique_parents())

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
        return locals()
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
        return locals()
    objfn_trans_kwargs = property(**objfn_trans_kwargs())

    def objfn_wt():
        doc = "The objfn_wt property."
        def fget(self):
            return self._objfn_wt
        def fset(self, value):
            self._objfn_wt = value
        def fdel(self):
            del self._objfn_wt
        return locals()
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
        return locals()
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
        return locals()
    ndset_trans_kwargs = property(**ndset_trans_kwargs())

    def ndset_wt():
        doc = "The ndset_wt property."
        def fget(self):
            return self._ndset_wt
        def fset(self, value):
            self._ndset_wt = value
        def fdel(self):
            del self._ndset_wt
        return locals()
    ndset_wt = property(**ndset_wt())

    def rng():
        doc = "The rng property."
        def fget(self):
            return self._rng
        def fset(self, value):
            if value is None:               # if None
                value = pybropt.core.random # use default random number generator
                return                      # exit function
            check_is_Generator(value, "rng")# check is numpy.Generator
            self._rng = value
        def fdel(self):
            del self._rng
        return locals()
    rng = property(**rng())

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
        mod : AdditiveLinearGenomicModel
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
        u           = mod.u                 # get regression coefficients

        if (chrgrp_stix is None) or (chrgrp_spix is None):
            raise RuntimeError("markers are not sorted by chromosome position")

        # get number of chromosomes
        nchr = len(chrgrp_stix)

        if self.nhaploblk < nchr:
            raise RuntimeError("number of haplotype blocks is less than the number of chromosomes")

        # calculate number of marker blocks to assign to each chromosome
        nblk = calc_nhaploblk_chrom(genpos, chrgrp_stix, chrgrp_spix)

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
    # FIXME: replace this function with select
    def pselect(self, t_cur, t_max, geno, bval, gmod, k = None, traitwt = None, **kwargs):
        """
        Select parents individuals for breeding.

        Parameters
        ----------
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        geno : dict
            A dict containing genotypic data for all breeding populations.
            Must have the following fields:
                Field | Type                         | Description
                ------+------------------------------+--------------------------
                cand  | PhasedGenotypeMatrix         | Parental candidate breeding population
                main  | PhasedGenotypeMatrix         | Main breeding population
                queue | List of PhasedGenotypeMatrix | Breeding populations on queue
                ""
        bval : dict
            A dict containing breeding value data.
            Must have the following fields:
                Field      | Type                        | Description
                -----------+-----------------------------+----------------------
                cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                main       | BreedingValueMatrix         | Main breeding population breeding values
                main_true  | BreedingValueMatrix         | Main breeding population true breeding values
        gmod : dict
            A dict containing genomic models.
            Must have the following fields:
                Field | Type                 | Description
                ------+----------------------+----------------------------------
                cand  | GenomicModel         | Parental candidate breeding population genomic model
                main  | GenomicModel         | Main breeding population genomic model
                true  | GenomicModel         | True genomic model for trait(s)
        k : int
        traitwt : numpy.ndarray
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing five objects: (pgvmat, sel, ncross, nprogeny, misc)
            pgvmat : PhasedGenotypeMatrix
                A PhasedGenotypeMatrix of parental candidates.
            sel : numpy.ndarray
                Array of indices specifying a cross pattern. Each index
                corresponds to an individual in 'pgvmat'.
            ncross : numpy.ndarray
                Number of crosses to perform per cross pattern.
            nprogeny : numpy.ndarray
                Number of progeny to generate per cross.
            misc : dict
                Miscellaneous output (user defined).
        """
        # get parameters
        if k is None:
            k = self.k_p
        if traitwt is None:
            traitwt = self.traitwt_p

        # construct pairs
        ntaxa = geno["cand"].ntaxa
        a = []
        for i in range(ntaxa):
            for j in range(i+1,ntaxa):
                a.append(i)
                a.append(j)
        a = numpy.int64(a)

        # get objective function
        objfn = self.pobjfn(
            t_cur = t_cur,
            t_max = t_max,
            geno = geno,
            bval = bval,
            gmod = gmod,
            traitwt = traitwt
        )

        gebv = objfn(a)                         # get all GEBVs
        if gebv.ndim == 1:                      # if there is one trait objective
            pairs = gebv.argsort()[::-1][:k]    # get indices of top k GEBVs
            self.rng.shuffle(pairs)             # shuffle indices
            sel = []
            for i in pairs:
                sel.append(a[2*i])      # female index
                sel.append(a[2*i+1])    # male index
            sel = numpy.int64(sel)
        elif gebv.ndim == 2:                    # TODO: ND-selection
            raise RuntimeError("non-dominated genomic selection not implemented")

        misc = {}

        return geno["cand"], sel, self.ncross, self.nprogeny, misc

    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, method = "single", nparent = None, ncross = None, nprogeny = None, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None, ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = None, **kwargs):
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
            A tuple containing five objects: (pgmat, sel, ncross, nprogeny, misc)
            pgmat : PhasedGenotypeMatrix
                A PhasedGenotypeMatrix of parental candidates.
            sel : numpy.ndarray
                Array of indices specifying a cross pattern. Each index
                corresponds to an individual in 'pgmat'.
            ncross : numpy.ndarray
                Number of crosses to perform per cross pattern.
            nprogeny : numpy.ndarray
                Number of progeny to generate per cross.
            misc : dict
                Miscellaneous output (user defined).
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

            # get all OHVs for each individual
            # (n,)
            ohv = [objfn([i]) for i in range(gmat.ntaxa)]

            # convert to numpy.ndarray
            ohv = numpy.array(ohv)

            # multiply the objectives by objfn_wt to transform to maximizing function
            # (n,) * scalar -> (n,)
            ohv = ohv * objfn_wt

            # get indices of top nparent GEBVs
            sel = ohv.argsort()[::-1][:nparent]

            # shuffle indices for random mating
            self.rng.shuffle(sel)

            # get GEBVs for reference
            misc = {"ohv" : ohv}

            return pgmat, sel, ncross, nprogeny, misc

        # multi-objective method: objfn_trans returns a multiple values for each
        # selection configuration
        elif method == "pareto":
            # get the pareto frontier
            frontier, sel_config, misc = self.pareto(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                nparent = nparent,
                objfn_trans = objfn_trans,
                objfn_trans_kwargs = objfn_trans_kwargs,
                objfn_wt = objfn_wt
            )

            # get scores for each of the points along the pareto frontier
            score = ndset_wt * ndset_trans(frontier, **ndset_trans_kwargs)

            # get index of maximum score
            ix = score.argmax()

            # add fields to misc
            misc["frontier"] = frontier
            misc["sel_config"] = sel_config

            return pgmat, sel_config[ix], ncross, nprogeny, misc
        else:
            raise ValueError("argument 'method' must be either 'single' or 'pareto'")

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, **kwargs):
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
        # get default parameters if any are None
        if trans is None:
            trans = self.objfn_trans
        if trans_kwargs is None:
            trans_kwargs = self.objfn_trans_kwargs

        mat = self._calc_hmat(pgmat, gpmod) # (m,n,b,t) get haplotype matrix
        xmap = self._calc_xmap()            # (s,p) get the cross map
        ploidy = pgmat.ploidy               # (scalar) get ploidy

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

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, **kwargs):
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
        # get default parameters if any are None
        if trans is None:
            trans = self.objfn_trans
        if trans_kwargs is None:
            trans_kwargs = self.objfn_trans_kwargs

        mat = self._calc_hmat(pgmat, gpmod) # (m,n,b,t) get haplotype matrix
        xmap = self._calc_xmap()            # (s,p) get the cross map
        ploidy = pgmat.ploidy               # (scalar) get ploidy

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

    # FIXME: add pareto function
    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, nparent = None, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None, **kwargs):
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
        # process inputs, apply defaults as needed.
        if nparent is None:
            nparent = self.nparent
        if objfn_trans is None:
            objfn_trans = self.objfn_trans
        if objfn_trans_kwargs is None:
            objfn_trans_kwargs = self.objfn_trans_kwargs
        if objfn_wt is None:
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

        # create optimization algorithm
        moalgo = NSGA2SetGeneticAlgorithm(
            rng = self.rng,
            **kwargs
        )

        # use multi-objective optimization to approximate Pareto front.
        frontier, sel_config, misc = moalgo.optimize(
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

                trans(numpy.ndarray, **kwargs):
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
        # (m,n,b,t)[:,(k,d),:,:] -> (m,k,d,b,t) # select k individuals
        # (m,k,d,b,t).max((0,2)) -> (k,b,t)     # find maximum haplotype across all parental phases
        # (k,b,t).sum((0,1)) -> (t,)            # add maximum haplotypes for k crosses and b blocks
        ohv = mat[:,sel,:,:].max((0,2)).sum((0,1))

        # multiply by ploidy
        # scalar * (t,) -> (t,)                 # multiply result by number of phases
        ohv *= ploidy

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

                trans(numpy.ndarray, **kwargs):
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
        # (m,n,b,t)[:,(j,k,d),:] -> (m,j,k,d,b,t)   # select k individuals
        # (m,j,k,d,b,t).max((0,3)) -> (j,k,b,t)     # find maximum haplotype across all parental phases
        # (j,k,b,t).sum((1,2)) -> (j,t)             # add maximum haplotypes for k crosses and b blocks
        ohv = mat[:,sel,:,:].max((0,3)).sum((1,2))

        # multiply by ploidy
        # scalar * (j,t) -> (j,t)                   # multiply result by number of phases
        ohv *= ploidy

        # apply transformations
        # (j,t) ---trans---> (?,?)
        if trans:
            ohv = trans(ohv, **kwargs)

        return ohv
