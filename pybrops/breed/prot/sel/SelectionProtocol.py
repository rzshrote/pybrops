"""
Module defining a general class for selection protocols.
"""

# list of all public objects in this module
__all__ = [
    "SelectionProtocol",
    "check_is_SelectionProtocol"
]

# imports
from abc import ABCMeta, abstractmethod
from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
from numpy.random import Generator, RandomState
import pandas

from pybrops.breed.prot.sel.cfg.SelectionConfiguration import SelectionConfiguration
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.prob.trans import trans_empty, trans_identity, trans_ndpt_to_vec_dist
from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState, check_is_ndarray, check_ndarray_dtype_is_integer
from pybrops.core.error.error_type_python import check_is_Callable, check_is_Integral, check_is_Real, check_is_dict
from pybrops.core.error.error_value_numpy import check_ndarray_all_gteq, check_ndarray_len_eq, check_ndarray_ndim
from pybrops.core.error.error_value_python import check_is_gt, check_is_gteq, check_is_neq
from pybrops.core.random.prng import global_prng
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class SelectionProtocol(metaclass=ABCMeta):
    """
    A semi-abstract class implementing several key properties common to most, 
    if not all, constrained selection protocols.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ncross: Integral,
            nparent: Integral,
            nmating: Union[Integral,numpy.ndarray],
            nprogeny: Union[Integral,numpy.ndarray],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            ndset_wt: Optional[Real] = None,
            ndset_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            rng: Optional[Union[Generator,RandomState]] = None, 
            soalgo: Optional[OptimizationAlgorithm] = None,
            moalgo: Optional[OptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the semi-abstract class SelectionProtocol.

        Parameters
        ----------
        ncross : Integral
            Number of cross configurations to consider.
        
        nparent : Integral
            Number of parents per cross configuration.
        
        nmating : Integral, numpy.ndarray
            Number of matings per configuration.

            If ``nmating`` is ``Integral``, then broadcast to a ``numpy.ndarray`` 
            of shape ``(ncross,)``.
            
            If ``nmating`` is ``numpy.ndarray``, then the array must be of type 
            ``Integral`` and of shape ``(ncross,)``.
        
        nprogeny : Integral, numpy.ndarray
            Number of progeny to derive from each mating event.

            If ``nprogeny`` is ``Integral``, then broadcast to a ``numpy.ndarray`` 
            of shape ``(ncross,)``.
            
            If ``nprogeny`` is ``numpy.ndarray``, then the array must be of type 
            ``Integral`` and of shape ``(ncross,)``.

        nobj : Integral
            Number of optimization objectives when constructing a 
            ``SelectionProblem``. This is equivalent to the vector length 
            returned by the ``obj_trans`` function. Must be ``Integral`` greater 
            than 0.
        
        obj_wt : numpy.ndarray, Real, None
            Objective function weights. Weights from this vector are applied 
            to objective function values via the Hadamard product. If values 
            are ``1.0`` or ``-1.0``, this can be used to specify minimizing 
            and maximizing objectives, respectively.

            If ``obj_wt`` is ``numpy.ndarray``, then the array must be of shape 
            ``(nobj,)``.

            If ``obj_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(nobj,)``.

            If ``obj_wt`` is ``None``, then the value ``1.0`` is broadcast to a 
            ``numpy.ndarray`` of shape ``(nobj,)``. This assumes that all 
            objectives are to be minimized.
        
        obj_trans : Callable, None
            A function which transforms values from a latent objective space to 
            the objective space. This transformation function must have the 
            following signature::

                def obj_trans(
                        decnvec: numpy.ndarray,
                        latentvec: numpy.ndarray, 
                        **kwargs: dict
                    ) -> numpy.ndarray:
                    # do stuff
                    return output
            
            Where:

            - ``decnvec`` is a ``numpy.ndarray`` containing the decision vector.
            - ``latentvec`` is a ``numpy.ndarray`` containing the latent space 
                objective function values which are to be transformed.
            - ``kwargs`` is a ``dict`` containing additional keyword arguments.

            If ``obj_trans`` is ``None``, then default to an identity objective 
            transformation function.

        obj_trans_kwargs : dict
            Keyword arguments for the latent space to objective space 
            transformation function. 

            If `obj_trans_kwargs`` is ``None``, then default to an empty 
            dictionary.

        nineqcv : Integral, None
            Number of inequality constraint violation functions. This is 
            equivalent to the vector length returned by the ``ineqcv_trans`` 
            function. Must be ``Integral`` greater than or equal to zero.

            If ``nineqcv`` is ``None``, then set to zero.

        ineqcv_wt : numpy.ndarray, None
            Inequality constraint violation function weights. Weights from this 
            vector are applied to inequality constraint violation function 
            values via the Hadamard product. If values are ``1.0`` or ``-1.0``, 
            this can be used to specify minimizing and maximizing constraints, 
            respectively.

            If ``ineqcv_wt`` is ``numpy.ndarray``, then the array must be of 
            shape ``(nineqcv,)``.

            If ``ineqcv_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(nineqcv,)``.

            If ``ineqcv_wt`` is ``None``, then the value ``1.0`` is broadcast 
            to a ``numpy.ndarray`` of shape ``(nineqcv,)``. This assumes that 
            all constraints are to be minimized.

        ineqcv_trans : Callable, None
            A function which transforms values from a latent objective space to 
            the inequality constraint violation space. This transformation 
            function must have the following signature::

                def ineqcv_trans(
                        decnvec: numpy.ndarray,
                        latentvec: numpy.ndarray, 
                        **kwargs: dict
                    ) -> numpy.ndarray:
                    # do stuff
                    return output
            
            Where:

            - ``decnvec`` is a ``numpy.ndarray`` containing the decision vector.
            - ``latentvec`` is a ``numpy.ndarray`` containing the latent space 
                objective function values which are to be transformed.
            - ``kwargs`` is a ``dict`` containing additional keyword arguments.

            If ``ineqcv_trans`` is ``None``, then default to a transformation 
            function returning an empty vector.
        
        ineqcv_trans_kwargs : dict, None
            Keyword arguments for the latent space to inequality constraint 
            violation transformation function.
        
            If `ineqcv_trans_kwargs`` is ``None``, then default to an empty 
            dictionary.

        neqcv : Integral, None
            Number of equality constraint violations. This is equivalent to the 
            vector length returned by the ``eqcv_trans`` function. Must be 
            ``Integral`` greater than or equal to zero.
        
            If ``neqcv`` is ``None``, then set to zero.

        eqcv_wt : numpy.ndarray, None
            Equality constraint violation function weights. Weights from this 
            vector are applied to equality constraint violation function 
            values via the Hadamard product. If values are ``1.0`` or ``-1.0``, 
            this can be used to specify minimizing and maximizing constraints, 
            respectively.

            If ``eqcv_wt`` is ``numpy.ndarray``, then the array must be of 
            shape ``(neqcv,)``.

            If ``eqcv_wt`` is ``Real``, then the value is broadcast to a 
            ``numpy.ndarray`` of shape ``(neqcv,)``.

            If ``eqcv_wt`` is ``None``, then the value ``1.0`` is broadcast 
            to a ``numpy.ndarray`` of shape ``(neqcv,)``. This assumes that 
            all constraints are to be minimized.

        eqcv_trans : Callable, None
            A function which transforms values from a latent objective space to 
            the equality constraint violation space. This transformation 
            function must have the following signature::

                def eqcv_trans(
                        decnvec: numpy.ndarray,
                        latentvec: numpy.ndarray, 
                        **kwargs: dict
                    ) -> numpy.ndarray:
                    # do stuff
                    return output
            
            Where:

            - ``decnvec`` is a ``numpy.ndarray`` containing the decision vector.
            - ``latentvec`` is a ``numpy.ndarray`` containing the latent space 
                objective function values which are to be transformed.
            - ``kwargs`` is a ``dict`` containing additional keyword arguments.

            If ``eqcv_trans`` is ``None``, then default to a transformation 
            function returning an empty vector.

        eqcv_trans_kwargs : dict, None
            Keyword arguments for the latent space to equality constraint 
            violation transformation function.

            If `eqcv_trans_kwargs`` is ``None``, then default to an empty 
            dictionary.

        ndset_wt : Real, None
            Nondominated set weight. The weight from this function is applied 
            to outputs from ``ndset_trans``. If values are ``1.0`` or ``-1.0``, 
            this can be used to specify minimizing and maximizing objectives, 
            respectively.

            If ``ndset_wt`` is ``None``, then it is set to the default value of ``1.0``.
            This assumes that the objective is to be minimized.

        ndset_trans : Callable, None
            A function which transforms values from the non-dominated set 
            objective space to the single-objective space. This transformation 
            function must have the following signature::

                def ndset_trans(
                        mat: numpy.ndarray, 
                        **kwargs: dict
                    ) -> numpy.ndarray:
                    # do stuff
                    return output
            
            Where:

            - ``mat`` is a ``numpy.ndarray`` containing a point coordinate array 
                of shape ``(npt, nobj)`` where ``npt`` is the number of points 
                and ``nobj`` is the number of objectives (dimensions). This 
                array contains input points for calculating the distance between 
                a point to the vector ``vec_wt``.
            - ``kwargs`` is a ``dict`` containing additional keyword arguments.

            If ``ndset_trans`` is ``None``, then default to a transformation 
            function calculating the distance between a weight vector and 
            provided points

        ndset_trans_kwargs : dict, None
            Nondominated set transformation function keyword arguments.

            If ``ndset_trans_kwargs`` is ``None``, then default to defaults for 
            the default ``ndset_trans`` function::

                ndset_trans_kwargs = {
                    "obj_wt": numpy.repeat(1.0, nobj),
                    "vec_wt": numpy.repeat(1.0, nobj)
                }

        rng : numpy.random.Generator, numpy.random.RandomState, None
            Random number source.

            If ``rng`` is ``None``, default to the global random number 
            generator.

        soalgo : OptimizationAlgorithm, None
            Single-objective optimization algorithm.

            If ``soalgo`` is ``None``, then use a default single-objective 
            optimization algorithm.

        moalgo : OptimizationAlgorithm, None
            Multi-objective opimization algorithm.

            If ``moalgo`` is ``None``, then use a default multi-objective 
            optimization algorithm.

        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        self.ncross = ncross
        self.nparent = nparent
        self.nmating = nmating
        self.nprogeny = nprogeny
        self.nobj = nobj
        self.obj_wt = obj_wt
        self.obj_trans = obj_trans
        self.obj_trans_kwargs = obj_trans_kwargs
        self.nineqcv = nineqcv
        self.ineqcv_wt = ineqcv_wt
        self.ineqcv_trans = ineqcv_trans
        self.ineqcv_trans_kwargs = ineqcv_trans_kwargs
        self.neqcv = neqcv
        self.eqcv_wt = eqcv_wt
        self.eqcv_trans = eqcv_trans
        self.eqcv_trans_kwargs = eqcv_trans_kwargs
        self.ndset_wt = ndset_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs
        self.rng = rng
        self.soalgo = soalgo
        self.moalgo = moalgo

    ############################ Object Properties #############################
    @property
    def nselindiv(self) -> Integral:
        """Number of selected individuals."""
        return self.ncross * self.nparent

    @property
    def ncross(self) -> Integral:
        """Number of cross configurations to consider."""
        return self._ncross
    @ncross.setter
    def ncross(self, value: Integral) -> None:
        """Set number of cross configurations."""
        check_is_Integral(value, "ncross")
        check_is_gt(value, "ncross", 0)
        self._ncross = value
    
    @property
    def nparent(self) -> Integral:
        """Number of parents per cross configuration."""
        return self._nparent
    @nparent.setter
    def nparent(self, value: Integral) -> None:
        """Set number of parents per cross configuration."""
        check_is_Integral(value, "nparent")
        check_is_gt(value, "nparent", 0)
        self._nparent = value

    @property
    def nmating(self) -> numpy.ndarray:
        """Number of matings per cross configuration."""
        return self._nmating
    @nmating.setter
    def nmating(self, value: Union[Integral,numpy.ndarray]) -> None:
        """Set number of matings per cross configuration."""
        if isinstance(value, Integral):
            check_is_gteq(value, "nmating", 0)
            value = numpy.repeat(value, self.ncross)
        check_is_ndarray(value, "nmating")
        check_ndarray_dtype_is_integer(value, "nmating")
        check_ndarray_all_gteq(value, "nmating", 0)
        self._nmating = value

    @property
    def nprogeny(self) -> numpy.ndarray:
        """Number of progeny to derive from each mating event."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: Union[Integral,numpy.ndarray]) -> None:
        """Set number of progeny to derive from each mating event."""
        if isinstance(value, Integral):
            check_is_gteq(value, "nprogeny", 0)
            value = numpy.repeat(value, self.ncross)
        check_is_ndarray(value, "nprogeny")
        check_ndarray_dtype_is_integer(value, "nprogeny")
        check_ndarray_all_gteq(value, "nprogeny", 0)
        self._nprogeny = value

    @property
    def nobj(self) -> Integral:
        """Number of optimization objectives."""
        return self._nobj
    @nobj.setter
    def nobj(self, value: Integral) -> None:
        """Set number of objectives."""
        check_is_Integral(value, "nobj")
        check_is_gt(value, "nobj", 0)     # cannot have 0 objectives
        self._nobj = value
    
    @property
    def obj_wt(self) -> numpy.ndarray:
        """Objective function weights."""
        return self._obj_wt
    @obj_wt.setter
    def obj_wt(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set objective function weights. If None, set to 1.0."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "obj_wt", 1)
            check_ndarray_len_eq(value, "obj_wt", self.nobj)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.nobj)
        elif value is None:
            value = numpy.repeat(1.0, self.nobj)
        else:
            raise TypeError("'obj_wt' must be of type numpy.ndarray, a real type, or None")
        self._obj_wt = value

    @property
    def obj_trans(self) -> Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to objective function values."""
        return self._obj_trans
    @obj_trans.setter
    def obj_trans(self, value: Union[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to objective space transformation function. If None, set to identity function."""
        if value is None:
            value = trans_identity
        check_is_Callable(value, "obj_trans")
        self._obj_trans = value
    
    @property
    def obj_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to objective space transformation function."""
        return self._obj_trans_kwargs
    @obj_trans_kwargs.setter
    def obj_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to objective space transformation function. If None, set to empty dict."""
        if value is None:
            value = {}
        check_is_dict(value, "obj_trans_kwargs")
        self._obj_trans_kwargs = value
    
    @property
    def nineqcv(self) -> Integral:
        """Number of inequality constraint violation functions."""
        return self._nineqcv
    @nineqcv.setter
    def nineqcv(self, value: Union[Integral,None]) -> None:
        """Set number of inequality constraint violation functions. If None, set to 0."""
        if value is None:
            value = 0
        check_is_Integral(value, "nineqcv")
        check_is_gteq(value, "nineqcv", 0)  # possible to have 0 inequality constraints
        self._nineqcv = value

    @property
    def ineqcv_wt(self) -> numpy.ndarray:
        """Inequality constraint violation function weights."""
        return self._ineqcv_wt
    @ineqcv_wt.setter
    def ineqcv_wt(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set inequality constraint violation function weights. If None, set to 1.0."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "ineqcv_wt", 1)
            check_ndarray_len_eq(value, "ineqcv_wt", self.nineqcv)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.nineqcv)
        elif value is None:
            value = numpy.repeat(1.0, self.nineqcv)
        else:
            raise TypeError("'ineqcv_wt' must be of type numpy.ndarray, a real type, or None")
        self._ineqcv_wt = value

    @property
    def ineqcv_trans(self) -> Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to inequality constraint violation values."""
        return self._ineqcv_trans
    @ineqcv_trans.setter
    def ineqcv_trans(self, value: Union[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to inequality constraint violation transformation function. If None, set to the empty function."""
        if value is None:
            value = trans_empty
        check_is_Callable(value, "ineqcv_trans")
        self._ineqcv_trans = value
    
    @property
    def ineqcv_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to inequality constraint violation transformation function."""
        return self._ineqcv_trans_kwargs
    @ineqcv_trans_kwargs.setter
    def ineqcv_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to inequality constraint violation transformation function. If None, set to empty dict."""
        if value is None:
            value = {}
        check_is_dict(value, "ineqcv_trans_kwargs")
        self._ineqcv_trans_kwargs = value
    
    @property
    def neqcv(self) -> Integral:
        """Number of equality constraint violations."""
        return self._neqcv
    @neqcv.setter
    def neqcv(self, value: Union[Integral,None]) -> None:
        """Set number of equality constraint violations. If None, set to 0."""
        if value is None:
            value = 0
        check_is_Integral(value, "neqcv")
        check_is_gteq(value, "neqcv", 0)    # possible to have 0 equality constraints
        self._neqcv = value
    
    @property
    def eqcv_wt(self) -> numpy.ndarray:
        """Equality constraint violation function weights."""
        return self._eqcv_wt
    @eqcv_wt.setter
    def eqcv_wt(self, value: Union[numpy.ndarray,Real,None]) -> None:
        """Set equality constraint violation function weights. If None, set to 1.0."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "eqcv_wt", 1)
            check_ndarray_len_eq(value, "eqcv_wt", self.neqcv)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.neqcv)
        elif value is None:
            value = numpy.repeat(1.0, self.neqcv)
        else:
            raise TypeError("'eqcv_wt' must be of type numpy.ndarray or a real type")
        self._eqcv_wt = value

    @property
    def eqcv_trans(self) -> Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]:
        """Function which transforms outputs from ``latentfn`` to equality constraint violation values."""
        return self._eqcv_trans
    @eqcv_trans.setter
    def eqcv_trans(self, value: Union[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to equality constraint violation transformation function. If None, set to the empty function."""
        if value is None:
            value = trans_empty
        check_is_Callable(value, "eqcv_trans")
        self._eqcv_trans = value 
    
    @property
    def eqcv_trans_kwargs(self) -> dict:
        """Keyword arguments for the latent space to equality constraint violation transformation function."""
        return self._eqcv_trans_kwargs
    @eqcv_trans_kwargs.setter
    def eqcv_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to equality constraint violation transformation function. If None, set to empty dict."""
        if value is None:
            value = {}
        check_is_dict(value, "eqcv_trans_kwargs")
        self._eqcv_trans_kwargs = value

    @property
    def ndset_wt(self) -> Real:
        """Nondominated set weights."""
        return self._ndset_wt
    @ndset_wt.setter
    def ndset_wt(self, value: Union[Real,None]) -> None:
        """Set nondominated set weights. If None, set to 1.0."""
        if value is None:
            value = 1.0
        check_is_Real(value, "ndset_wt")
        check_is_neq(value, "ndset_wt", 0.0)
        self._ndset_wt = value

    @property
    def ndset_trans(self) -> Callable[[numpy.ndarray,dict],numpy.ndarray]:
        """Nondominated set transformation function."""
        return self._ndset_trans
    @ndset_trans.setter
    def ndset_trans(self, value: Union[Callable[[numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set nondominated set transformation function. If None, set to closest point to vector function."""
        if value is None:
            value = trans_ndpt_to_vec_dist
        check_is_Callable(value, "ndset_trans")
        self._ndset_trans = value

    @property
    def ndset_trans_kwargs(self) -> dict:
        """Nondominated set transformation function keyword arguments."""
        return self._ndset_trans_kwargs
    @ndset_trans_kwargs.setter
    def ndset_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set nondominated set transformation function keyword arguments. If None, set to dict with one vectors."""
        if value is None:                           # if given None
            value = {                               # set default to empty dict
                "obj_wt": numpy.repeat(1.0, self.nobj),
                "vec_wt": numpy.repeat(1.0, self.nobj)
            }
        check_is_dict(value, "ndset_trans_kwargs")  # check is dict
        self._ndset_trans_kwargs = value

    @property
    def rng(self) -> Union[Generator,RandomState]:
        """rng."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[Generator,RandomState,None]) -> None:
        """Set rng."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng")
        self._rng = value

    @property
    @abstractmethod
    def soalgo(self) -> OptimizationAlgorithm:
        """Single-objective optimization algorithm."""
        raise NotImplementedError("property is abstract")
    @soalgo.setter
    @abstractmethod
    def soalgo(self, value: Union[OptimizationAlgorithm,None]) -> None:
        """Set single-objective optimization algorithm."""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def moalgo(self) -> OptimizationAlgorithm:
        """Multi-objective opimization algorithm."""
        raise NotImplementedError("property is abstract")
    @moalgo.setter
    @abstractmethod
    def moalgo(self, value: Union[OptimizationAlgorithm,None]) -> None:
        """Set multi-objective opimization algorithm."""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    @abstractmethod
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: pandas.DataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            **kwargs: dict
        ) -> SelectionProblem:
        """
        Create an optimization problem definition using provided inputs.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : pandas.DataFrame
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
        out : SelectionProblem
            An optimization problem definition.
        """
        raise NotImplementedError("method is abstract")

    ################ Single Objective Solve ################
    @abstractmethod
    def sosolve(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: pandas.DataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict], 
            **kwargs: dict
        ) -> SelectionSolution:
        """
        Solve the selection problem using a single-objective optimization algorithm.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : pandas.DataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Solution
            A single-objective solution to the posed selection problem.
        """
        raise NotImplementedError("method is abstract")

    ################ Multi Objective Solve #################
    @abstractmethod
    def mosolve(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: pandas.DataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict], 
            **kwargs: dict
        ) -> SelectionSolution:
        """
        Solve the selection problem using a multi-objective optimization algorithm.
        This calculates a Pareto frontier for the objectives.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : pandas.DataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Solution
            A multi-objective solution to the posed selection problem.
        """
        raise NotImplementedError("method is abstract")

    ################# Selection Functions ##################
    @abstractmethod
    def select(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: pandas.DataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict], 
            **kwargs: dict
        ) -> SelectionConfiguration:
        """
        Select a single selection configuration of individuals for breeding.
        
        If there is a single objective, then optimize using a single-objective 
        optimization algorithm and create a selection configuration from the
        identified solution.

        If there is are multiple objectives, then optimize using a multi-objective
        optimization algorithm, apply a transformation on the non-dominated set
        of solutions to find the single best solution among the set, and create
        a selection configuration from the identified solution.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : pandas.DataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : SelectionConfiguration
            A selection configuration object, requiring all necessary information to mate individuals.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_SelectionProtocol(v: object, vname: str) -> None:
    """
    Check if object is of type ConstrainedSelectionProtocol, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SelectionProtocol):
        raise TypeError("'{0}' must be of type ConstrainedSelectionProtocol.".format(vname))
