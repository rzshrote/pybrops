import numpy
import types

import pybropt.core.random

from pybropt.breed.prot.sel.SelectionProtocol import SelectionProtocol

from pybropt.core.error import check_is_bool
from pybropt.core.error import check_is_int
from pybropt.core.error import cond_check_is_Generator
from pybropt.core.error import cond_check_is_dict
from pybropt.core.error import cond_check_is_callable

class RandomSelection(SelectionProtocol):
    """Perform random parent selection"""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, nparent, ncross, nprogeny, replace = False, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = 1.0, ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = 1.0, rng = None, **kwargs):
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
            If ``rng`` is ``None``, use ``pybropt.core.random`` module
            (NOT THREAD SAFE!).
        """
        super(RandomSelection, self).__init__(**kwargs)

        # check inputs
        check_is_int(nparent, "nparent")
        check_is_int(ncross, "ncross")
        check_is_int(nprogeny, "nprogeny")
        check_is_bool(replace, "replace")
        cond_check_is_callable(objfn_trans, "objfn_trans")
        cond_check_is_dict(objfn_trans_kwargs, "objfn_trans_kwargs")
        # TODO: check objfn_wt
        cond_check_is_callable(ndset_trans, "ndset_trans")
        cond_check_is_dict(ndset_trans_kwargs, "ndset_trans_kwargs")
        # TODO: check ndset_wt
        cond_check_is_Generator(rng, "rng")

        # assign variables
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.replace = replace
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = {} if objfn_trans_kwargs is None else objfn_trans_kwargs
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = {} if ndset_trans_kwargs is None else ndset_trans_kwargs
        self.ndset_wt = ndset_wt
        self.rng = pybropt.core.random if rng is None else rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, nparent = None, ncross = None, nprogeny = None, replace = None, **kwargs):
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

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, **kwargs):
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

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, **kwargs):
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

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs):
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
        rng : numpy.Generator
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
