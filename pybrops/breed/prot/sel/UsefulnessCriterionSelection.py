"""
Module implementing selection protocols for Usefulness Criterion selection.
"""

import numbers
import types
from typing import Callable, Optional, Type, Union
import numpy
import scipy.stats

from pybrops.algo.opt.NSGA2SetGeneticAlgorithm import NSGA2SetGeneticAlgorithm
from pybrops.algo.opt.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.algo.opt.SteepestAscentSetHillClimber import SteepestAscentSetHillClimber
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.core.error import check_is_bool
from pybrops.core.error import check_is_callable
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_gt
from pybrops.core.error import check_is_str
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_type_python import check_is_Number
from pybrops.core.error.error_value_python import check_is_gteq, check_is_lt
from pybrops.core.random.prng import global_prng
from pybrops.core.util.arrayix import triuix
from pybrops.core.util.arrayix import triudix
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory, check_is_GeneticVarianceMatrixFactory
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class UsefulnessCriterionSelection(SelectionProtocol):
    """
    Class implementing selection protocols for usefulness criterion selection.

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
            nself: int,
            upper_percentile: numbers.Number,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool = True, 
            method: str = "single",
            objfn_trans: Optional[Callable] = None, 
            objfn_trans_kwargs: Optional[dict] = None, 
            objfn_wt: Union[numpy.ndarray,numbers.Number] = 1.0,
            ndset_trans: Optional[Callable] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            ndset_wt: Union[numpy.ndarray,numbers.Number] = 1.0,
            rng: Union[numpy.random.Generator,numpy.random.RandomState] = None, 
            soalgo: Optional[OptimizationAlgorithm] = None, 
            moalgo: Optional[OptimizationAlgorithm] = None,
            **kwargs : dict
        ):
        """
        Constructor for Usefulness Criterion Selection (UCS).

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
            | ``"single"`` | UCS is transformed to a single objective and      |
            |              | optimization is done on the transformed function. |
            |              | This is done using the ``trans`` function         |
            |              | provided::                                        |
            |              |                                                   |
            |              |    optimize : objfn_trans(MOGS)                   |
            +--------------+---------------------------------------------------+
            | ``"pareto"`` | UCS is transformed by a transformation function,  |
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
        super(UsefulnessCriterionSelection, self).__init__(**kwargs)

        # error checks and assignments (ORDER DEPENDENT!!!)
        self.nconfig = nconfig
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.nself = nself
        self.upper_percentile = upper_percentile
        self.vmatfcty = vmatfcty
        self.gmapfn = gmapfn
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
        """Description for property nparent."""
        return self._nparent
    @nparent.setter
    def nparent(self, value: int) -> None:
        """Set data for property nparent."""
        check_is_int(value, "nparent")      # must be int
        check_is_gt(value, "nparent", 0)    # int must be >0
        self._nparent = value
    @nparent.deleter
    def nparent(self) -> None:
        """Delete data for property nparent."""
        del self._nparent
    
    @property
    def ncross(self) -> int:
        """Description for property ncross."""
        return self._ncross
    @ncross.setter
    def ncross(self, value: int) -> None:
        """Set data for property ncross."""
        check_is_int(value, "ncross")       # must be int
        check_is_gt(value, "ncross", 0)     # int must be >0
        self._ncross = value
    @ncross.deleter
    def ncross(self) -> None:
        """Delete data for property ncross."""
        del self._ncross
    
    @property
    def nprogeny(self) -> int:
        """Description for property nprogeny."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: int) -> None:
        """Set data for property nprogeny."""
        check_is_int(value, "nprogeny")     # must be int
        check_is_gt(value, "nprogeny", 0)   # int must be >0
        self._nprogeny = value
    @nprogeny.deleter
    def nprogeny(self) -> None:
        """Delete data for property nprogeny."""
        del self._nprogeny
    
    @property
    def nself(self) -> int:
        """Description for property nself."""
        return self._nself
    @nself.setter
    def nself(self, value: int) -> None:
        """Set data for property nself."""
        check_is_int(value, "nself")     # must be int
        check_is_gteq(value, "nself", 0)   # int must be >=0
        self._nself = value
    @nself.deleter
    def nself(self) -> None:
        """Delete data for property nself."""
        del self._nself
    
    @property
    def unique_parents(self) -> bool:
        """Description for property unique_parents."""
        return self._unique_parents
    @unique_parents.setter
    def unique_parents(self, value: bool) -> None:
        """Set data for property unique_parents."""
        check_is_bool(value, "unique_parents")
        self._unique_parents = value
    @unique_parents.deleter
    def unique_parents(self) -> None:
        """Delete data for property unique_parents."""
        del self._unique_parents
    
    @property
    def upper_percentile(self) -> numbers.Number:
        """Description for property upper_percentile."""
        return self._upper_percentile
    @upper_percentile.setter
    def upper_percentile(self, value: numbers.Number) -> None:
        """Set data for property upper_percentile."""
        check_is_Number(value, "upper_percentile")  # must be a number
        check_is_gt(value, "upper_percentile", 0)   # number must be >0
        check_is_lt(value, "upper_percentile", 1)   # number must be <1
        self._upper_percentile = value
        # set the selection intensity
        self._selection_intensity = scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - self._upper_percentile)) / self._upper_percentile
    @upper_percentile.deleter
    def upper_percentile(self) -> None:
        """Delete data for property upper_percentile."""
        del self._upper_percentile
    
    @property
    def selection_intensity(self) -> numbers.Number:
        """Description for property selection_intensity."""
        return self._selection_intensity
    @selection_intensity.setter
    def selection_intensity(self, value: numbers.Number) -> None:
        """Set data for property selection_intensity."""
        error_readonly("selection_intensity")
    @selection_intensity.deleter
    def selection_intensity(self) -> None:
        """Delete data for property selection_intensity."""
        error_readonly("selection_intensity")

    @property
    def vmatfcty(self) -> GeneticVarianceMatrixFactory:
        """Description for property vmatfcty."""
        return self._vmatfcty
    @vmatfcty.setter
    def vmatfcty(self, value: GeneticVarianceMatrixFactory) -> None:
        """Set data for property vmatfcty."""
        check_is_GeneticVarianceMatrixFactory(value, "vmatfcty")
        self._vmatfcty = value
    @vmatfcty.deleter
    def vmatfcty(self) -> None:
        """Delete data for property vmatfcty."""
        del self._vmatfcty
    
    @property
    def gmapfn(self) -> GeneticMapFunction:
        """Description for property gmapfn."""
        return self._gmapfn
    @gmapfn.setter
    def gmapfn(self, value: GeneticMapFunction) -> None:
        """Set data for property gmapfn."""
        self._gmapfn = value
    @gmapfn.deleter
    def gmapfn(self) -> None:
        """Delete data for property gmapfn."""
        del self._gmapfn

    @property
    def method(self) -> str:
        """Description for property method."""
        return self._method
    @method.setter
    def method(self, value: str) -> None:
        """Set data for property method."""
        check_is_str(value, "method")       # must be string
        value = value.lower()               # convert to lowercase
        options = ("single", "pareto")      # method options
        if value not in options:            # if not method supported
            # raise ValueError
            raise ValueError("Unsupported method. Options are: " + ", ".join(map(str, options)))
        self._method = value
    @method.deleter
    def method(self) -> None:
        """Delete data for property method."""
        del self._method
    
    @property
    def objfn_trans(self) -> Callable:
        """Description for property objfn_trans."""
        return self._objfn_trans
    @objfn_trans.setter
    def objfn_trans(self, value: Union[Callable,None]) -> None:
        """Set data for property objfn_trans."""
        if value is not None:                       # if given object
            check_is_callable(value, "objfn_trans") # must be callable
        self._objfn_trans = value
    @objfn_trans.deleter
    def objfn_trans(self) -> None:
        """Delete data for property objfn_trans."""
        del self._objfn_trans
    
    @property
    def objfn_trans_kwargs(self) -> dict:
        """Description for property objfn_trans_kwargs."""
        return self._objfn_trans_kwargs
    @objfn_trans_kwargs.setter
    def objfn_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set data for property objfn_trans_kwargs."""
        if value is None:                           # if given None
            value = {}                              # set default to empty dict
        check_is_dict(value, "objfn_trans_kwargs")  # check is dict
        self._objfn_trans_kwargs = value
    @objfn_trans_kwargs.deleter
    def objfn_trans_kwargs(self) -> None:
        """Delete data for property objfn_trans_kwargs."""
        del self._objfn_trans_kwargs
    
    @property
    def objfn_wt(self) -> numpy.ndarray:
        """Description for property objfn_wt."""
        return self._objfn_wt
    @objfn_wt.setter
    def objfn_wt(self, value: numpy.ndarray) -> None:
        """Set data for property objfn_wt."""
        self._objfn_wt = value
    @objfn_wt.deleter
    def objfn_wt(self) -> None:
        """Delete data for property objfn_wt."""
        del self._objfn_wt
    
    @property
    def ndset_trans(self) -> Union[Callable,None]:
        """Description for property ndset_trans."""
        return self._ndset_trans
    @ndset_trans.setter
    def ndset_trans(self, value: Union[Callable,None]) -> None:
        """Set data for property ndset_trans."""
        if value is not None:                       # if given object
            check_is_callable(value, "ndset_trans") # must be callable
        self._ndset_trans = value
    @ndset_trans.deleter
    def ndset_trans(self) -> None:
        """Delete data for property ndset_trans."""
        del self._ndset_trans
    
    @property
    def ndset_trans_kwargs(self) -> dict:
        """Description for property ndset_trans_kwargs."""
        return self._ndset_trans_kwargs
    @ndset_trans_kwargs.setter
    def ndset_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set data for property ndset_trans_kwargs."""
        if value is None:                           # if given None
            value = {}                              # set default to empty dict
        check_is_dict(value, "ndset_trans_kwargs")  # check is dict
        self._ndset_trans_kwargs = value
    @ndset_trans_kwargs.deleter
    def ndset_trans_kwargs(self) -> None:
        """Delete data for property ndset_trans_kwargs."""
        del self._ndset_trans_kwargs
    
    @property
    def ndset_wt(self) -> numpy.ndarray:
        """Description for property ndset_wt."""
        return self._ndset_wt
    @ndset_wt.setter
    def ndset_wt(self, value: numpy.ndarray) -> None:
        """Set data for property ndset_wt."""
        self._ndset_wt = value
    @ndset_wt.deleter
    def ndset_wt(self) -> None:
        """Delete data for property ndset_wt."""
        del self._ndset_wt
    
    @property
    def rng(self) -> Union[numpy.random.Generator,numpy.random.RandomState]:
        """Description for property rng."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[numpy.random.Generator,numpy.random.RandomState]) -> None:
        """Set data for property rng."""
        # if None, use default random number generator
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng")# check is numpy.Generator
        self._rng = value
    @rng.deleter
    def rng(self) -> None:
        """Delete data for property rng."""
        del self._rng
    
    @property
    def soalgo(self) -> OptimizationAlgorithm:
        """Description for property soalgo."""
        return self._soalgo
    @soalgo.setter
    def soalgo(self, value: Union[OptimizationAlgorithm,None]) -> None:
        """Set data for property soalgo."""
        if value is None:
            value = SteepestAscentSetHillClimber(
                rng = self.rng  # PRNG source
            )
        self._soalgo = value
    @soalgo.deleter
    def soalgo(self) -> None:
        """Delete data for property soalgo."""
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
        self._moalgo = value
    @moalgo.deleter
    def moalgo(self) -> None:
        """Delete multi-objective opimization algorithm."""
        del self._moalgo
    
    ############################################################################
    ########################## Private Object Methods ##########################
    ############################################################################
    def _calc_xmap(self, ntaxa: int) -> numpy.ndarray:
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

    def _calc_uc(self, pgmat: PhasedGenotypeMatrix, gmod: GenomicModel, xmap: numpy.ndarray):
        # calculate breeding values
        bvmat_obj = gmod.gebv(pgmat)
        
        # calculate variance matrix
        vmat_obj = self.vmatfcty.from_gmod(
            gmod = gmod, 
            pgmat = pgmat, 
            ncross = self.ncross, 
            nprogeny = self.nprogeny, 
            nself = self.nself,
            gmapfn = self.gmapfn
        )

        # get expected genome contributions
        # (1,p)
        epgc = numpy.ndarray(vmat_obj.epgc)
        
        # extract decentered and descaled
        bvmat = bvmat_obj.descale()     # (n,t)

        # extract variances
        vmat = vmat_obj.mat             # (n,...,n,t)

        # allocate memory for usefulness criterion
        uc = numpy.empty((len(xmap),bvmat_obj.ntrait), dtype = float)

        # for each cross configuration
        for i,cconfig in enumerate(xmap):
            # take dot product with expected genome contributions to get progeny mean
            # (p,) . (p,t) -> (t,)
            pmean = epgc.dot(bvmat[cconfig,:])

            # extract the variance at (i,...,i,:)
            # (t,)
            pvar = vmat[tuple(cconfig) + (slice(None),)]

            # calculate the usefulness criterion
            # (t,) + (t,) -> (t,)
            uc[i,:] = pmean + self.selection_intensity * numpy.sqrt(pvar)
        
        return uc

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
        # get selection parameters
        nconfig = self.nconfig
        nparent = self.nparent
        ncross = self.ncross
        nprogeny = self.nprogeny
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
                **kwargs
            )

            # calculate xmap
            xmap = self._calc_xmap(pgmat.ntaxa)

            # get all UCSs for each configuration
            # (s,)
            ucs = [objfn([i]) for i in range(len(xmap))]

            # convert to numpy.ndarray
            ucs = numpy.array(ucs)

            # multiply the objectives by objfn_wt to transform to maximizing function
            # (n,) * scalar -> (n,)
            ucs = ucs * objfn_wt

            # get indices of top nconfig UCSs
            sel = ucs.argsort()[::-1][:nconfig]

            # shuffle indices for random mating
            self.rng.shuffle(sel)

            # convert 'sel' to parent selections (ordered)
            # (kd,)
            sel = xmap[sel,:].flatten()

            # get GEBVs for reference
            misc = {"UCS" : ucs}

            # add optimization details to miscellaneous output
            if miscout is not None:     # if miscout was provided
                miscout.update(misc)    # add dict to dict

            return pgmat, sel, ncross, nprogeny

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
                **kwargs
            )

            # get scores for each of the points along the pareto frontier
            score = ndset_wt * ndset_trans(frontier, **ndset_trans_kwargs)

            # get index of maximum score
            ix = score.argmax()

            # add fields to misc
            misc["frontier"] = frontier
            misc["sel_config"] = sel_config

            return pgmat, sel_config[ix], ncross, nprogeny, misc

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
        xmap = self._calc_xmap(pgmat.ntaxa)     # (s,p) get the cross map
        uc = self._calc_uc(pgmat, gpmod, xmap)  # (s,t) get usefulness criterion matrix
        trans = self.objfn_trans                # get transformation function
        trans_kwargs = self.objfn_trans_kwargs  # get transformation function keyword arguments

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,         # byte code pointer
            self.objfn_static.__globals__,      # global variables
            None,                               # new name for the function
            (uc, trans, trans_kwargs),          # default values for arguments
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
        xmap = self._calc_xmap(pgmat.ntaxa)     # (s,p) get the cross map
        uc = self._calc_uc(pgmat, gpmod, xmap)  # (s,t) get usefulness criterion matrix
        trans = self.objfn_trans                # get transformation function
        trans_kwargs = self.objfn_trans_kwargs  # get transformation function keyword arguments

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,     # byte code pointer
            self.objfn_vec_static.__globals__,  # global variables
            None,                               # new name for the function
            (uc, trans, trans_kwargs),          # default values for arguments
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
        objfn_wt = self.objfn_wt
        moalgo = self.moalgo

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

        # TODO: fix sel_config output format: currently in xmap indices format.
        return frontier, sel_config

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, uc, trans, kwargs):
        """
        Score a population of individuals based on Usefulness Criterion
        Selection (UCS). Scoring for UCS is defined as the sum of maximum
        haploid breeding values obtainable from a population.

        UCS selects the ``q`` individuals with the largest UCSs.

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
        UCS : numpy.ndarray
            A UCS matrix of shape ``(t,)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.
        """
        # get the cross configurations
        # (s,t)[(k,),:] -> (k,t)
        # (k,t).sum(0) -> (t,)
        ucs = uc[sel,:].sum(0)

        # apply transformations
        # (t,) ---trans---> (?,)
        if trans:
            ucs = trans(ucs, **kwargs)

        return ucs

    @staticmethod
    def objfn_vec_static(sel, uc, trans, kwargs):
        """
        Score a population of individuals based on Usefulness Criterion
        Selection (UCS). Scoring for UCS is defined as the sum of maximum
        haploid breeding values obtainable from a population.

        UCS selects the ``q`` individuals with the largest UCSs.

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
        UCS : numpy.ndarray
            A UCS matrix of shape ``(t,)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.
        """
        # get the cross configurations
        # (s,t)[(j,k),:] -> (j,k,t)
        # (j,k,t).sum(0) -> (j,t)
        ucs = uc[sel,:].sum(1)

        # apply transformations
        # (j,t) ---trans---> (j,?)
        if trans:
            ucs = trans(ucs, **kwargs)

        return ucs
