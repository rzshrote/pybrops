"""
Module implementing selection protocols for multi-objective genomic selection.
"""

import numpy
import math
import types
from typing import Union

import pybrops.core.random
from pybrops.algo.opt.NSGA2SetGeneticAlgorithm import NSGA2SetGeneticAlgorithm
from pybrops.algo.opt.SteepestAscentSetHillClimber import SteepestAscentSetHillClimber
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.core.error import check_isinstance
from pybrops.core.error import check_is_bool
from pybrops.core.error import check_is_callable
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_gt
from pybrops.core.error import check_is_str
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.util.pareto import is_pareto_efficient
from pybrops.core.random.prng import global_prng
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel

class PopulationAlleleUnavailabilitySelection(SelectionProtocol):
    """
    Class implementing selection protocols for multi-objective genomic selection.

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
        target: Union[str,numpy.ndarray] = "positive", 
        weight: Union[str,numpy.ndarray] = "magnitude", 
        method: str = "single",
        objfn_trans = None, 
        objfn_trans_kwargs = None, 
        objfn_wt = -1.0,
        ndset_trans = None, 
        ndset_trans_kwargs = None, 
        ndset_wt = -1.0,
        soalgo = None, 
        moalgo = None,
        rng = global_prng, 
        **kwargs
        ):
        """
        Constructor for PopulationAlleleUnavailabilitySelection class.

        Parameters
        ----------
        nparent : int
            Number of parents to select.
        ncross : int
            Number of crosses per configuration.
        nprogeny : int
            Number of progeny to derive from each cross.
        target : str or numpy.ndarray
            If target is a string, check value and follow these rules:

            +-------------------+----------------------------------------------+
            | Value             | Description                                  |
            +===================+==============================================+
            | ``"positive"``    | Select alleles with the most positive        |
            |                   | effect.                                      |
            +-------------------+----------------------------------------------+
            | ``"negative"``    | Select alleles with the most negate effect.  |
            +-------------------+----------------------------------------------+
            | ``"stabilizing"`` | Set target allele frequency to ``0.5``.      |
            +-------------------+----------------------------------------------+
            | ``numpy.ndarray`` | Use frequency values in ``target`` as is.    |
            +-------------------+----------------------------------------------+
        weight : str or numpy.ndarray
            If weight is a string, check value and follow these rules:

            +-----------------+------------------------------------------------+
            | Value           | Description                                    |
            +=================+================================================+
            | ``"magnitude"`` | Assign weights using the magnitudes of         |
            |                 | regression coefficients.                       |
            +-----------------+------------------------------------------------+
            | ``"equal"``     | Assign weights equally.                        |
            +-----------------+------------------------------------------------+
        method : str
            Method of selecting parents.

            +--------------+---------------------------------------------------+
            | Method       | Description                                       |
            +==============+===================================================+
            | ``"single"`` | MOGS is transformed to a single objective and     |
            |              | optimization is done on the transformed function. |
            |              | This is done using the ``trans`` function         |
            |              | provided::                                        |
            |              |                                                   |
            |              |    optimize : objfn_trans(MOGS)                   |
            +--------------+---------------------------------------------------+
            | ``"pareto"`` | MOGS is transformed by a transformation function, |
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
        objfn_trans : function or callable
            Function to transform the MOGS function. If method = "single", this
            function must return a scalar. If method = "pareto", this function
            must return a numpy.ndarray.

            Function definition::

                objfn_trans(obj, **kwargs):
                    Parameters
                        obj : scalar, numpy.ndarray
                            Objective scalar or vector to be transformed
                        kwargs : dict
                            Additional keyword arguments
                    Returns
                        out : scalar, numpy.ndarray
                            Transformed objective scalar or vector.
        objfn_trans_kwargs : dict
            Dictionary of keyword arguments to be passed to 'objfn_trans'.
        objfn_wt : float, numpy.ndarray
            Weight applied to transformed objective function. Indicates whether
            a function is maximizing or minimizing:

            - ``1.0`` for maximizing function.
            - ``-1.0`` for minimizing function.
        ndset_trans : numpy.ndarray
            Function to transform nondominated points along the Pareto frontier
            into a single score for each point.

            Function definition::

                ndset_trans(ndset, **kwargs):
                    Parameters
                        ndset : numpy.ndarray
                            Array of shape (j,o) containing nondominated points.
                            Where 'j' is the number of nondominated points and
                            'o' is the number of objectives.
                        kwargs : dict
                            Additional keyword arguments.
                    Returns
                        out : numpy.ndarray
                            Array of shape (j,) containing transformed Pareto
                            frontier points.
        ndset_trans_kwargs : dict
            Dictionary of keyword arguments to be passed to 'ndset_trans'.
        ndset_wt : float
            Weight applied to transformed nondominated points along Pareto
            frontier. Indicates whether a function is maximizing or minimizing.
                1.0 for maximizing function.
                -1.0 for minimizing function.
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
        rng : numpy.random.Generator or None
            A random number generator source. Used for optimization algorithms.
            If 'rng' is None, use pybrops.core.random module (NOT THREAD SAFE!).
        """
        super(PopulationAlleleUnavailabilitySelection, self).__init__(**kwargs)

        # error checks and assignments (ORDER DEPENDENT!!!)
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.target = target
        self.weight = weight
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

    def target():
        doc = "The target property."
        def fget(self):
            return self._target
        def fset(self, value):
            check_isinstance(value, "target", (str, numpy.ndarray))
            if isinstance(value, str):
                value = value.lower()               # convert to lowercase
                options = (                         # target options
                    'positive',
                    'negative',
                    'stabilizing'
                )
                if value not in options:            # if target not supported
                    raise ValueError(               # raise ValueError
                        "Unsupported 'target'. Options are: " +
                        ", ".join(map(str, options))
                    )
            self._target = value
        def fdel(self):
            del self._target
        return locals()
    target = property(**target())

    def weight():
        doc = "The weight property."
        def fget(self):
            return self._weight
        def fset(self, value):
            check_isinstance(value, "weight", (str, numpy.ndarray))
            if isinstance(value, str):
                value = value.lower()               # convert to lowercase
                options = ('magnitude', 'equal')    # weight options
                if value not in options:            # if weight not supported
                    raise ValueError(               # raise ValueError
                        "Unsupported 'weight'. Options are: " +
                        ", ".join(map(str, options))
                    )
            self._weight = value
        def fdel(self):
            del self._weight
        return locals()
    weight = property(**weight())

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
        return locals()
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
        return locals()
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
        return locals()
    moalgo = property(**moalgo())

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
        return locals()
    rng = property(**rng())

    ############################################################################
    ########################## Private Object Methods ##########################
    ############################################################################
    def _calc_mkrwt(self, gpmod: AdditiveLinearGenomicModel):
        if isinstance(self.weight, str):
            if self.weight == "magnitude":
                return numpy.absolute(gpmod.u_a)
            elif self.weight == "equal":
                return 1.0
                # return numpy.full(gpmod.u_a.shape, 1.0, dtype = "float64")
            else:
                raise ValueError("string value for 'weight' not recognized")
        elif isinstance(self.weight, numpy.ndarray):
            return self.weight
        else:
            raise TypeError("variable 'weight' must be a string or numpy.ndarray")

    def _calc_tfreq(self, gpmod: AdditiveLinearGenomicModel):
        if isinstance(self.target, str):
            if self.target == "positive":
                return numpy.float64(gpmod.u_a >= 0.0)   # positive alleles are desired
            elif self.target == "negative":
                return numpy.float64(gpmod.u_a <= 0.0)   # negative alleles are desired
            elif self.target == "stabilizing":
                return 0.5                          # both alleles desired
                # return numpy.full(coeff.shape, 0.5, dtype = 'float64')
            else:
                raise ValueError("string value for 'target' not recognized")
        elif isinstance(self.target, numpy.ndarray):
            return self.target
        else:
            raise TypeError("variable 'target' must be a string or numpy.ndarray")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs):
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
        method = self.method
        objfn_trans = self.objfn_trans
        objfn_trans_kwargs = self.objfn_trans_kwargs
        objfn_wt = self.objfn_wt
        ndset_trans = self.ndset_trans
        ndset_trans_kwargs = self.ndset_trans_kwargs
        ndset_wt = self.ndset_wt

        # single-objective method: objfn_trans returns a single value for each
        # selection configuration
        if method == "single":
            # get number of taxa
            ntaxa = pgmat.ntaxa

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

            # optimize using hill-climber algorithm
            opt = self.soalgo.optimize(
                objfn,                          # objective function
                k = nparent,                    # number of parents to select
                sspace = numpy.arange(ntaxa),   # parental indices
                objfn_wt = objfn_wt,            # maximizing function
                **kwargs
            )

            # get best solution
            sel = opt["soln"]

            # add optimization details to miscellaneous output
            if miscout is not None:     # if miscout was provided
                miscout.update(opt)     # add dict to dict

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
                **kwargs
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

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
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
        gpmod : AdditiveLinearGenomicModel
            Linear genomic prediction model.

        Returns
        -------
        outfn : function
            A selection objective function for the specified problem.
        """
        # get selection parameters
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        # calculate default function parameters
        mat = gmat.mat                      # (n,p) get genotype matrix
        ploidy = gmat.ploidy                # (scalar) get number of phases
        mkrwt = self._calc_mkrwt(gpmod)     # (p,t) get marker weights
        tfreq = self._calc_tfreq(gpmod)     # (p,t) get target allele frequencies

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,                         # byte code pointer
            self.objfn_static.__globals__,                      # global variables
            None,                                               # new name for the function
            (mat, ploidy, tfreq, mkrwt, trans, trans_kwargs),   # default values for arguments
            self.objfn_static.__closure__                       # closure byte code pointer
        )

        return outfn

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
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
        gpmod : AdditiveLinearGenomicModel
            Linear genomic prediction model.

        Returns
        -------
        outfn : function
            A vectorized selection objective function for the specified problem.
        """
        # get selection parameters
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs
        weight = self.weight
        target = self.target

        # calculate default function parameters
        mat = gmat.mat                      # (n,p) get genotype matrix
        ploidy = gmat.ploidy                # (scalar) get number of phases
        u = gpmod.u_a                       # (p,t) get regression coefficients
        mkrwt = self._calc_mkrwt(weight, u) # (p,t) get marker weights
        tfreq = self._calc_tfreq(target, u) # (p,t) get target allele frequencies

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,                     # byte code pointer
            self.objfn_vec_static.__globals__,                  # global variables
            None,                                               # new name for the function
            (mat, ploidy, tfreq, mkrwt, trans, trans_kwargs),   # default values for arguments
            self.objfn_vec_static.__closure__                   # closure byte code pointer
        )

        return outfn

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs):
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
    def objfn_static(sel, mat, ploidy, tfreq, mkrwt, trans, kwargs):
        """
        Population Allele Unavailability (PAU): :math:`f^{\\textup{PAU}}(\\textbf{x})`

        .. math::

            f^{\\textup{PAU}}(\\textbf{x}) = \\textbf{w} \\cdot \\textbf{u}

        Given a genotype matrix ``mat`` and a selection indices vector
        :math:`\\textbf{x} =` ``sel``, calculate the selection allele frequency.
        From the selection allele frequencies and the target allele frequencies
        ``tfreq``, determine if the target frequencies can be attained after
        unlimited generations of selection. If the target allele frequency at a
        locus cannot be attained, score locus as ``1``, otherwise score as
        ``0``. Store this into a binary score vector :math:`\\textbf{u}`.
        Take the dot product between the binary score vector and the marker
        weight vector :math:`\\textbf{w} =` ``mkrwt`` to calculate
        :math:`f^{\\textup{PAU}}(\\textbf{x})` and return the result.

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
            A genotype matrix of shape ``(n,p)`` representing only biallelic
            loci. One of the two alleles at a locus is coded using a ``1``. The
            other allele is coded as a ``0``. ``mat`` holds the counts of the
            allele coded by ``1``.

            Where:

            - ``n`` is the number of individuals.
            - ``p`` is the number of markers.

            Example::

                # matrix of shape (n = 3, p = 4)
                mat = numpy.array([[0,2,1,0],
                                   [2,2,1,1],
                                   [0,1,0,2]])
        ploidy : int
            Number of phases that the genotype matrix ``mat`` represents.
        tfreq : floating, numpy.ndarray
            A target allele frequency matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Example::

                tfreq = numpy.array([0.2, 0.6, 0.7, 0.5])
        mkrwt : numpy.ndarray
            A marker weight coefficients matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Remarks:

            - All values in ``mkrwt`` must be non-negative.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single numpy.ndarray argument.
            - Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        mogs : numpy.ndarray
            A MOGS score matrix of shape ``(t + t,)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.

            Matrix element ordering for un-transformed MOGS score matrix:

            - The first set of ``t`` elements in the ``mogs`` output correspond
              to the ``t`` PAU outputs for each trait.
            - The second set of ``t`` elements in the ``mogs`` output correspond
              to the ``t`` PAFD outputs for each trait.
        """
        # if no selection, select all
        if sel is None:
            sel = slice(None)

        # generate a view of the genotype matrix that only contains 'sel' rows.
        # (n,p)[(k,),:] -> (k,p)
        sgeno = mat[sel,:]

        # calculate reciprocal number of phases
        # number of phases * number of individuals in 'sgeno'
        rphase = 1.0 / (ploidy * sgeno.shape[1])

        # calculate population frequencies; add axis for correct broadcast
        # (k,p).sum(0) -> (p,)
        # (p,) * scalar -> (p,)
        # (p,None) -> (p,1)
        # Remark: we need (p,1) for broadcasting with (p,t) arrays
        pfreq = (sgeno.sum(0) * rphase)[:,None]

        # calculate some inequalities for use multiple times
        pfreq_lteq_0 = (pfreq <= 0.0)   # is population freq <= 0.0
        pfreq_gteq_1 = (pfreq >= 1.0)   # is population freq >= 1.0

        # calculate allele unavailability
        # (p,t)
        allele_unavail = numpy.where(
            tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
            pfreq_lteq_0,           # then set True if sel has allele freq == 0
            numpy.where(            # else
                tfreq > 0.0,        # if 0.0 < target freq < 1.0
                numpy.logical_or(   # then set True if pop freq is outside (0.0,1.0)
                    pfreq_lteq_0,   # mask for whether population freq <= 0.0
                    pfreq_gteq_1    # mask for whether population freq >= 1.0
                ),
                pfreq_gteq_1        # else set True if pop freq is >= 1.0
            )
        )

        # compute PAU(x)
        # (p,t) * (p,t) -> (p,t)
        # (p,t).sum[0] -> (t,)
        pau = (mkrwt * allele_unavail).sum(0)

        # apply transformations
        if trans is not None:
            pau = trans(pau, **kwargs)

        return pau

    @staticmethod
    def objfn_vec_static(sel, mat, ploidy, tfreq, mkrwt, trans, kwargs):
        """
        A vectorized multi-objective genomic selection objective function.

        - The goal is to minimize all objectives for this function.
        - This is a bare bones function. Minimal error checking is done.

        Objectives: :math:`F(\\textbf{x})`

        .. math::

            F(\\textbf{x}) = {[f^{\\textup{PAU}}(\\textbf{x}), f^{\\textup{PAFD}}(\\textbf{x})]}'

        Population Allele Unavailability (PAU): :math:`f^{\\textup{PAU}}(\\textbf{x})`

        .. math::

            f^{\\textup{PAU}}(\\textbf{x}) = \\textbf{w} \\cdot \\textbf{u}

        Given a genotype matrix ``mat`` and a selection indices vector
        :math:`\\textbf{x} =` ``sel``, calculate the selection allele frequency.
        From the selection allele frequencies and the target allele frequencies
        ``tfreq``, determine if the target frequencies can be attained after
        unlimited generations of selection. If the target allele frequency at a
        locus cannot be attained, score locus as ``1``, otherwise score as
        ``0``. Store this into a binary score vector :math:`\\textbf{u}`.
        Take the dot product between the binary score vector and the marker
        weight vector :math:`\\textbf{w} =` ``mkrwt`` to calculate
        :math:`f^{\\textup{PAU}}(\\textbf{x})` and return the result.

        Population Allele Frequency Distance (PAFD): :math:`f^{\\textup{PAFD}}(\\textbf{x})`

        .. math::
            f^{\\textup{PAFD}}(\\textbf{x}) = \\textbf{w} \\cdot \\left | \\textbf{p}_{x} - \\textbf{p}_{t} \\right |

        Given a genotype matrix ``mat`` and a selection indices vector
        :math:`\\textbf{x} =` ``sel``, calculate the selection allele frequency
        :math:`\\textbf{p}_{x}`. From the selection allele frequencies and the
        target allele frequencies :math:`\\textbf{p}_{t} =` ``tfreq``,
        calculate the absolute value of the difference between the two vectors.
        Finally, take the dot product between the difference vector and the marker
        weight vector :math:`\\textbf{w} =` ``mkrwt`` to calculate
        :math:`f^{\\textup{PAFD}}(\\textbf{x})` and return the result.

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
        mat : numpy.ndarray
            A genotype matrix of shape ``(n,p)`` representing only biallelic
            loci. One of the two alleles at a locus is coded using a ``1``. The
            other allele is coded as a ``0``. ``mat`` holds the counts of the
            allele coded by ``1``.

            Where:

            - ``n`` is the number of individuals.
            - ``p`` is the number of markers.

            Example::

                # matrix of shape (n = 3, p = 4)
                mat = numpy.array([[0,2,1,0],
                                   [2,2,1,1],
                                   [0,1,0,2]])
        ploidy : int
            Number of phases that the genotype matrix ``mat`` represents.
        tfreq : floating, numpy.ndarray
            A target allele frequency matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Example::

                tfreq = numpy.array([0.2, 0.6, 0.7, 0.5])
        mkrwt : numpy.ndarray
            A marker weight coefficients matrix of shape ``(p,t)``.

            Where:

            - ``p`` is the number of markers.
            - ``t`` is the number of traits.

            Remarks:

            - All values in ``mkrwt`` must be non-negative.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        mogs : numpy.ndarray
            A MOGS score matrix of shape ``(j,t + t)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``t`` is the number of traits.

            Matrix element ordering for un-transformed MOGS score matrix:

            - The first set of ``t`` elements in the ``mogs`` output correspond
              to the ``t`` PAU outputs for each trait.
            - The second set of ``t`` elements in the ``mogs`` output correspond
              to the ``t`` PAFD outputs for each trait.
        """
        # generate a view of the genotype matrix that only contains 'sel' rows.
        # (n,p)[(j,k),:] -> (j,k,p)
        sgeno = mat[sel,:]

        # calculate reciprocal number of phases
        # ploidy * number of individuals in 'sgeno'
        rphase = 1.0 / (ploidy * sgeno.shape[2])

        # calculate population frequencies; add axis for correct broadcast
        # (j,k,p).sum(1) -> (j,p)
        # (j,p) * scalar -> (j,p)
        # (j,p)[:,None] -> (j,p,1)
        # Remark: we need (j,p,1) for broadcasting with (p,t) arrays
        pfreq = (sgeno.sum(1) * rphase)[:,None]

        # calculate some inequalities for use multiple times
        pfreq_lteq_0 = (pfreq <= 0.0)   # (j,p,1) is population freq <= 0.0
        pfreq_gteq_1 = (pfreq >= 1.0)   # (j,p,1) is population freq >= 1.0

        # calculate allele unavailability
        # (j,p,t)
        allele_unavail = numpy.where(
            tfreq >= 1.0,           # (p,t) if target freq >= 1.0 (should always be 1.0)
            pfreq_lteq_0,           # (j,p,1) then set True if sel has allele freq == 0
            numpy.where(            # (j,p,t) else
                tfreq > 0.0,        # (p,t) if 0.0 < target freq < 1.0
                numpy.logical_or(   # (j,p,1) then set True if pop freq is outside (0.0,1.0)
                    pfreq_lteq_0,   # (j,p,1) mask for whether population freq <= 0.0
                    pfreq_gteq_1    # (j,p,1) mask for whether population freq >= 1.0
                ),
                pfreq_gteq_1        # (j,p,1) else set True if pop freq is >= 1.0
            )
        )

        # compute f_PAU(x)
        # (p,t) * (j,p,t) -> (j,p,t)
        # (j,p,t).sum[1] -> (j,t)
        pau = (mkrwt * allele_unavail).sum(1)

        # apply transformations
        if trans is not None:
            pau = trans(pau, **kwargs)

        return pau
