import numpy
import types
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol

import pybrops.core.random
from pybrops.core.error import cond_check_is_Generator
from pybrops.core.error import check_is_int
from pybrops.core.error import cond_check_is_callable
from pybrops.core.error import cond_check_is_dict
from pybrops.algo.opt.NSGA2SetGeneticAlgorithm import NSGA2SetGeneticAlgorithm

class ConventionalPhenotypicSelection(SelectionProtocol):
    """docstring for ConventionalPhenotypicSelection."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, nparent, ncross, nprogeny, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = 1.0, ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = 1.0, rng = None, **kwargs):
        """
        Constructor for Conventional Phenotypic Selection (CPS)

        Parameters
        ----------
        nparent : int
        ncross : int
        nprogeny : int
        objfn_trans : function, callable, None
        objfn_trans_kwargs : dict, None
        objfn_wt : float, numpy.ndarray
        ndset_trans : function, callable, None
        ndset_trans_kwargs : dict, None
        ndset_wt : float
        rng : numpy.Generator
        """
        super(ConventionalPhenotypicSelection, self).__init__(**kwargs)

        # error checks
        check_is_int(nparent, "nparent")
        check_is_int(ncross, "ncross")
        check_is_int(nprogeny, "nprogeny")
        cond_check_is_callable(objfn_trans, "objfn_trans")
        cond_check_is_dict(objfn_trans_kwargs, "objfn_trans_kwargs")
        # TODO: check objfn_wt
        cond_check_is_callable(ndset_trans, "ndset_trans")
        cond_check_is_dict(ndset_trans_kwargs, "ndset_trans_kwargs")
        # TODO: check ndset_wt
        cond_check_is_Generator(rng, "rng")

        # variable assignment
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = {} if objfn_trans_kwargs is None else objfn_trans_kwargs
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = {} if ndset_trans_kwargs is None else ndset_trans_kwargs
        self.ndset_wt = ndset_wt
        self.rng = pybrops.core.random if rng is None else rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, method = "single", nparent = None, ncross = None, nprogeny = None, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None, ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = None, **kwargs):
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

            # get all EBVs for each individual
            # (n,)
            ebv = [objfn([i]) for i in range(bvmat.ntaxa)]

            # convert to numpy.ndarray
            ebv = numpy.array(ebv)

            # multiply the objectives by objfn_wt to transform to maximizing function
            # (n,) * scalar -> (n,)
            ebv = ebv * objfn_wt

            # get indices of top nparent GEBVs
            sel = ebv.argsort()[::-1][:nparent]

            # shuffle indices for random mating
            self.rng.shuffle(sel)

            # get GEBVs for reference
            if miscout is not None:
                miscout["ebv"] = ebv

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

        # raise error since method is not supported
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
            Not used by this function.
        ptdf : PhenotypeDataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Used by this function. Input breeding value matrix.
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

        # get pointers to raw numpy.ndarray matrices
        mat = bvmat.mat     # (n,t) get breeding value matrix

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,     # byte code pointer
            self.objfn_static.__globals__,  # global variables
            None,                           # new name for the function
            (mat, trans, trans_kwargs),     # default values for arguments
            self.objfn_static.__closure__   # closure byte code pointer
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
            Not used by this function.
        ptdf : PhenotypeDataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Used by this function. Input breeding value matrix.
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

        # get pointers to raw numpy.ndarray matrices
        mat = bvmat.mat     # (n,t) get breeding value matrix

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,     # byte code pointer
            self.objfn_vec_static.__globals__,  # global variables
            None,                               # new name for the function
            (mat, trans, trans_kwargs),         # default values for arguments
            self.objfn_vec_static.__closure__   # closure byte code pointer
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
        ntaxa = bvmat.ntaxa

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

        frontier, sel_config, misc = moalgo.optimize(
            objfn = objfn,                  # objective function
            k = nparent,                    # vector length to optimize (sspace^k)
            sspace = numpy.arange(ntaxa),   # search space options
            objfn_wt = objfn_wt             # weights to apply to each objective
        )

        if miscout is not None:
            for k,i in misc.items():
                miscout[k] = i

        return frontier, sel_config

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, mat, trans, kwargs):
        """
        Score a selection configuration based on its breeding values
        (Conventional Phenotype Selection; CPS).

        CPS selects the ``q`` individuals with the largest EBVs.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            If ``sel`` is ``None``, use all individuals.
        mat : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        cps : numpy.ndarray
            A EBV matrix of shape ``(t,)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.
        """
        # if sel is None, slice all individuals
        if sel is None:
            sel = slice(None)

        # CPS calculation explanation
        # Step 1: (n,t) -> (k,t)        # select individuals
        # Step 2: (k,t).sum(0) -> (t,)  # sum across all individuals
        cps = mat[sel,:].sum(0)

        # apply transformations
        # (t,) ---trans---> (?,)
        if trans:
            cps = trans(cps, **kwargs)

        return cps

    @staticmethod
    def objfn_vec_static(sel, mat, trans, kwargs):
        """
        Score a selection configuration based on its breeding values
        (Conventional Phenotype Selection; CPS).

        CPS selects the ``q`` individuals with the largest EBVs.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape ``(j,k)``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            If ``sel`` is ``None``, score each individual separately: ``(n,1)``
        mat : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        cps : numpy.ndarray
            A GEBV matrix of shape ``(j,t)`` if ``trans`` is None.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``t`` is the number of traits.
        """
        # if sel is None, slice all individuals
        if sel is None:
            n = mat.shape[0]
            sel = numpy.arange(n).reshape(n,1)

        # CPS calculation explanation
        # Step 1: (n,t)[(j,k),:] -> (j,k,t) # select individuals
        # Step 2: (j,k,t).sum(1) -> (j,t)   # sum across all individuals
        cps = mat[sel,:].sum(1)

        # apply transformations
        # (j,t) ---trans---> (?,?)
        if trans:
            cps = trans(cps, **kwargs)

        return cps