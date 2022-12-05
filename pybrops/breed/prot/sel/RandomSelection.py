"""
Module implementing selection protocols for random selection.
"""

import numpy
import types

import pybrops.core.random
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.core.error import check_is_bool
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_callable
from pybrops.core.random.prng import global_prng

class RandomSelection(SelectionProtocol):
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

    def replace():
        doc = "The replace property."
        def fget(self):
            """Get value for replace."""
            return self._replace
        def fset(self, value):
            """Set value for replace."""
            check_is_bool(value, "replace")
            self._replace = value
        def fdel(self):
            """Delete value for replace."""
            del self._replace
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    replace = property(**replace())

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
