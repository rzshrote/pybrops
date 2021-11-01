import numpy
import types

from . import SelectionProtocol

import pybropt.core.random
from pybropt.core.error import check_is_int
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import cond_check_is_ndarray
from pybropt.core.error import cond_check_is_Generator
from pybropt.core.error import cond_check_is_callable
from pybropt.core.error import cond_check_is_dict
from pybropt.algo.opt import NSGA2SetGeneticAlgorithm

class WeightedGenomicSelection(SelectionProtocol):
    """docstring for WeightedGenomicSelection."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, nparent, ncross, nprogeny, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = 1.0, ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = 1.0, rng = None, **kwargs):
        """
        Constructor for WeightedGenomicSelection class.

        Parameters
        ----------
        """
        super(WeightedGenomicSelection, self).__init__(**kwargs)

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
        self.rng = pybropt.core.random if rng is None else rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, method = "single", nparent = None, ncross = None, nprogeny = None, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None, ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = None, **kwargs):
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
        method : str
            Options: "single", "pareto"
        nparent : int, None
            Number of parents. If None, use default.
        ncross : int, None
            Number of crosses per configuration. If None, use default.
        nprogeny : int
            Number of progeny per cross. If None, use default.
        **kwargs
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
            misc = {"wgebv" : wgebv}

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
        Return an objective function for the provided datasets.
        """
        # get default parameters if any are None
        if trans is None:
            trans = self.objfn_trans
        if trans_kwargs is None:
            trans_kwargs = self.objfn_trans_kwargs

        # get pointers to raw numpy.ndarray matrices
        mat = gmat.mat  # (n,p) get genotype matrix
        u = gpmod.u     # (p,t) get regression coefficients

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
        """
        # get default parameters if any are None
        if trans is None:
            trans = self.objfn_trans
        if trans_kwargs is None:
            trans_kwargs = self.objfn_trans_kwargs

        # get pointers to raw numpy.ndarray matrices
        mat = gmat.mat      # (n,p) get genotype matrix
        u = gpmod.u   # (p,t) get regression coefficients

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

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, nparent = None, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None, **kwargs):
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
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing three objects (frontier, sel_config, misc)
            Elements
            --------
            frontier : numpy.ndarray
                Array of shape (q,v) containing Pareto frontier points.
                Where:
                    'q' is the number of points in the frontier.
                    'v' is the number of objectives for the frontier.
            sel_config : numpy.ndarray
                Array of shape (q,k) containing parent selection decisions for
                each corresponding point in the Pareto frontier.
                Where:
                    'q' is the number of points in the frontier.
                    'k' is the number of search space decision variables.
            misc : dict
                A dictionary of miscellaneous output. (User specified)
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

        frontier, sel_config, misc = moalgo.optimize(
            objfn = objfn,                  # objective function
            k = nparent,                    # vector length to optimize (sspace^k)
            sspace = numpy.arange(ntaxa),   # search space options
            objfn_wt = objfn_wt             # weights to apply to each objective
        )

        return frontier, sel_config, misc

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, mat, u, uwt, trans, kwargs):
        """
        Score a population of individuals based on Weighted Genomic Selection
        (WGS). Scoring for WGS is defined as the sum of weighted Genomic
        Estimated Breeding Values (wGEBV) for a population.

        WGS selects the 'q' individuals with the largest GEBVs.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape (k,)
            Where:
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, use all individuals.
        mat : numpy.ndarray
            A int8 binary genotype matrix of shape (n, p).
            Where:
                'n' is the number of individuals.
                'p' is the number of markers.
        u : numpy.ndarray
            A trait prediction coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        uwt : numpy.ndarray
            Multiplicative marker weights matrix to apply to the trait
            prediction coefficients provided of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            Trait prediction coefficients (u) are transformed as follows:
                u_new = u ⊙ uwt (Hadamard product)
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:
                Must accept a single numpy.ndarray argument.
                Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to 'trans' function.

        Returns
        -------
        wgs : numpy.ndarray
            A GEBV matrix of shape (k, t) if objwt is None.
            A GEBV matrix of shape (k,) if objwt shape is (t,)
            Where:
                'k' is the number of individuals selected.
                't' is the number of traits.
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
            A selection indices matrix of shape (j,k)
            Where:
                'j' is the number of selection configurations.
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, score each individual separately: (n,1)
        mat : numpy.ndarray
            A genotype matrix of shape (n, p).
            Where:
                'n' is the number of individuals.
                'p' is the number of markers.
        u : numpy.ndarray
            A trait prediction coefficients matrix of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
        uwt : numpy.ndarray
            Multiplicative marker weights matrix to apply to the trait
            prediction coefficients provided of shape (p, t).
            Where:
                'p' is the number of markers.
                't' is the number of traits.
            Trait prediction coefficients (u) are transformed as follows:
                u_new = u ⊙ uwt (Hadamard product)
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:
                Must accept a single numpy.ndarray argument.
                Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to 'trans' function.

        Returns
        -------
        cgs : numpy.ndarray
            A GEBV matrix of shape (j,t) if 'trans' is None.
            Otherwise, of shape specified by 'trans'.
            Where:
                'j' is the number of selection configurations.
                't' is the number of traits.
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
