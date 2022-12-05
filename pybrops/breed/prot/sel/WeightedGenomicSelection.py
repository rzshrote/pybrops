"""
Module implementing selection protocols for weighted genomic selection.
"""

import numpy
import types

import pybrops.core.random
from pybrops.algo.opt.NSGA2SetGeneticAlgorithm import NSGA2SetGeneticAlgorithm
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_ndarray
from pybrops.core.error import check_is_callable
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.random.prng import global_prng

class WeightedGenomicSelection(SelectionProtocol):
    """
    Class implementing selection protocols for weighted genomic selection.

    # TODO: add formulae for methodology.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, nparent, ncross, nprogeny, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = 1.0, ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = 1.0, rng = None, **kwargs):
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
    def nparent():
        doc = "The nparent property."
        def fget(self):
            """Get value for nparent."""
            return self._nparent
        def fset(self, value):
            """Set value for nparent."""
            check_is_int(value, "nparent")
            self._nparent = value
        def fdel(self):
            """Delete value for nparent."""
            del self._nparent
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    nparent = property(**nparent())

    def ncross():
        doc = "The ncross property."
        def fget(self):
            """Get value for ncross."""
            return self._ncross
        def fset(self, value):
            """Set value for ncross."""
            check_is_int(value, "ncross")
            self._ncross = value
        def fdel(self):
            """Delete value for ncross."""
            del self._ncross
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    ncross = property(**ncross())

    def nprogeny():
        doc = "The nprogeny property."
        def fget(self):
            """Get value for nprogeny."""
            return self._nprogeny
        def fset(self, value):
            """Set value for nprogeny."""
            check_is_int(value, "nprogeny")
            self._nprogeny = value
        def fdel(self):
            """Delete value for nprogeny."""
            del self._nprogeny
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    nprogeny = property(**nprogeny())

    def objfn_trans():
        doc = "The objfn_trans property."
        def fget(self):
            """Get value for objfn_trans."""
            return self._objfn_trans
        def fset(self, value):
            """Set value for objfn_trans."""
            if value is not None:
                check_is_callable(value, "objfn_trans")
            self._objfn_trans = value
        def fdel(self):
            """Delete value for objfn_trans."""
            del self._objfn_trans
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    objfn_trans = property(**objfn_trans())

    def objfn_trans_kwargs():
        doc = "The objfn_trans_kwargs property."
        def fget(self):
            """Get value for objfn_trans_kwargs."""
            return self._objfn_trans_kwargs
        def fset(self, value):
            """Set value for objfn_trans_kwargs."""
            if value is not None:
                check_is_dict(value, "objfn_trans_kwargs")
            else:
                value = {}
            self._objfn_trans_kwargs = value
        def fdel(self):
            """Delete value for objfn_trans_kwargs."""
            del self._objfn_trans_kwargs
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    objfn_trans_kwargs = property(**objfn_trans_kwargs())

    # TODO: finish error checks
    def objfn_wt():
        doc = "The objfn_wt property."
        def fget(self):
            """Get value for objfn_wt."""
            return self._objfn_wt
        def fset(self, value):
            """Set value for objfn_wt."""
            self._objfn_wt = value
        def fdel(self):
            """Delete value for objfn_wt."""
            del self._objfn_wt
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    objfn_wt = property(**objfn_wt())

    def ndset_trans():
        doc = "The ndset_trans property."
        def fget(self):
            """Get value for ndset_trans."""
            return self._ndset_trans
        def fset(self, value):
            """Set value for ndset_trans."""
            if value is not None:
                check_is_callable(value, "ndset_trans")
            self._ndset_trans = value
        def fdel(self):
            """Delete value for ndset_trans."""
            del self._ndset_trans
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    ndset_trans = property(**ndset_trans())

    def ndset_trans_kwargs():
        doc = "The ndset_trans_kwargs property."
        def fget(self):
            """Get value for ndset_trans_kwargs."""
            return self._ndset_trans_kwargs
        def fset(self, value):
            """Set value for ndset_trans_kwargs."""
            if value is not None:
                check_is_dict(value, "ndset_trans_kwargs")
            else:
                value = {}
            self._ndset_trans_kwargs = value
        def fdel(self):
            """Delete value for ndset_trans_kwargs."""
            del self._ndset_trans_kwargs
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    ndset_trans_kwargs = property(**ndset_trans_kwargs())

    # TODO: finish error checks
    def ndset_wt():
        doc = "The ndset_wt property."
        def fget(self):
            """Get value for ndset_wt."""
            return self._ndset_wt
        def fset(self, value):
            """Set value for ndset_wt."""
            self._ndset_wt = value
        def fdel(self):
            """Delete value for ndset_wt."""
            del self._ndset_wt
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    ndset_wt = property(**ndset_wt())

    def rng():
        doc = "The rng property."
        def fget(self):
            """Get value for rng."""
            return self._rng
        def fset(self, value):
            """Set value for rng."""
            if value is not None:
                check_is_Generator_or_RandomState(value, "rng")
            else:
                value = global_prng
            self._rng = value
        def fdel(self):
            """Delete value for rng."""
            del self._rng
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    rng = property(**rng())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, method = "single", nparent = None, ncross = None, nprogeny = None, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None, ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = None, **kwargs):
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

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, **kwargs):
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

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, **kwargs):
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
