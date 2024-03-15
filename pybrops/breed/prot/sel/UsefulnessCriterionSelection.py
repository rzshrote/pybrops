"""
Module implementing selection protocols for Usefulness Criterion (UC) selection.
"""

__all__ = [
    "UsefulnessCriterionSelectionMixin",
    "UsefulnessCriterionBinarySelection",
    "UsefulnessCriterionIntegerSelection",
    "UsefulnessCriterionRealSelection",
    "UsefulnessCriterionSubsetSelection",
]

from abc import ABCMeta
from numbers import Integral, Real
from typing import Optional, Union
from typing import Callable

import numpy
from numpy.random import Generator, RandomState
import pandas
import scipy.stats

from pybrops.breed.prot.sel.BinaryMateSelectionProtocol import BinaryMateSelectionProtocol
from pybrops.breed.prot.sel.IntegerMateSelectionProtocol import IntegerMateSelectionProtocol
from pybrops.breed.prot.sel.RealMateSelectionProtocol import RealMateSelectionProtocol
from pybrops.breed.prot.sel.SubsetMateSelectionProtocol import SubsetMateSelectionProtocol
from pybrops.breed.prot.sel.prob.BinaryMateSelectionProblem import BinaryMateSelectionProblem
from pybrops.breed.prot.sel.prob.IntegerMateSelectionProblem import IntegerMateSelectionProblem
from pybrops.breed.prot.sel.prob.RealMateSelectionProblem import RealMateSelectionProblem
from pybrops.breed.prot.sel.prob.SubsetMateSelectionProblem import SubsetMateSelectionProblem
from pybrops.breed.prot.sel.prob.UsefulnessCriterionSelectionProblem import UsefulnessCriterionBinaryMateSelectionProblem
from pybrops.breed.prot.sel.prob.UsefulnessCriterionSelectionProblem import UsefulnessCriterionIntegerMateSelectionProblem
from pybrops.breed.prot.sel.prob.UsefulnessCriterionSelectionProblem import UsefulnessCriterionRealMateSelectionProblem
from pybrops.breed.prot.sel.prob.UsefulnessCriterionSelectionProblem import UsefulnessCriterionSubsetMateSelectionProblem
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory, check_is_GeneticVarianceMatrixFactory
from pybrops.opt.algo.BinaryOptimizationAlgorithm import BinaryOptimizationAlgorithm
from pybrops.opt.algo.IntegerOptimizationAlgorithm import IntegerOptimizationAlgorithm
from pybrops.core.error.error_type_python import check_is_Integral, check_is_Real, check_is_bool
from pybrops.core.error.error_value_python import check_is_gt, check_is_gteq, check_is_in_interval_exclusive
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.opt.algo.RealOptimizationAlgorithm import RealOptimizationAlgorithm
from pybrops.opt.algo.SubsetOptimizationAlgorithm import SubsetOptimizationAlgorithm
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction, check_is_GeneticMapFunction
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix

class UsefulnessCriterionSelectionMixin(
        metaclass = ABCMeta,
    ):
    """
    Semi-abstract class for Usefulness Criterion (UC) Selection with constraints.
    """

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################
    @property
    def ntrait(self) -> Integral:
        """Number of traits to expect from breeding value matrix inputs."""
        return self._ntrait
    @ntrait.setter
    def ntrait(self, value: Integral) -> None:
        """Set number of traits to expect."""
        check_is_Integral(value, "ntrait")
        check_is_gt(value, "ntrait", 0)
        self._ntrait = value

    @property
    def nself(self) -> Integral:
        """Number of selfing generations post-cross before progenies are derived. Used in variance calculations."""
        return self._nself
    @nself.setter
    def nself(self, value: Integral) -> None:
        """Set the number of selfing generations post-cross before progenies are derived."""
        check_is_Integral(value, "nself")     # must be int
        check_is_gteq(value, "nself", 0)   # int must be >=0
        self._nself = value

    @property
    def upper_percentile(self) -> Real:
        """Upper percentile of progenies to select."""
        return self._upper_percentile
    @upper_percentile.setter
    def upper_percentile(self, value: Real) -> None:
        """Set upper percentile of progenies to select."""
        check_is_Real(value, "upper_percentile")  # must be a number
        check_is_in_interval_exclusive(value, "upper_percentile", 0.0, 1.0)
        self._upper_percentile = value
    
    @property
    def selection_intensity(self) -> Real:
        """Get selection intensity."""
        return scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - self._upper_percentile)) / self._upper_percentile

    @property
    def vmatfcty(self) -> GeneticVarianceMatrixFactory:
        """Variance matrix factory to calculate progeny genetic variances."""
        return self._vmatfcty
    @vmatfcty.setter
    def vmatfcty(self, value: GeneticVarianceMatrixFactory) -> None:
        """Set variance matrix factory to calculate progeny genetic variances."""
        check_is_GeneticVarianceMatrixFactory(value, "vmatfcty")
        self._vmatfcty = value
    
    @property
    def gmapfn(self) -> GeneticMapFunction:
        """Genetic map function used in the calculation of progeny genetic variance matrices."""
        return self._gmapfn
    @gmapfn.setter
    def gmapfn(self, value: GeneticMapFunction) -> None:
        """Set genetic map function used in the calculation of progeny genetic variance matrices."""
        check_is_GeneticMapFunction(value, "gmapfn")
        self._gmapfn = value

    @property
    def unique_parents(self) -> bool:
        """Whether parents should be unique."""
        return self._unique_parents
    @unique_parents.setter
    def unique_parents(self, value: bool) -> None:
        """Set whether parents should be unique."""
        check_is_bool(value, "unique_parents")
        self._unique_parents = value

class UsefulnessCriterionBinarySelection(UsefulnessCriterionSelectionMixin,BinaryMateSelectionProtocol):
    """
    Class defining Usefulness Criterion (UC) Selection for a binary search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool,
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
            soalgo: Optional[BinaryOptimizationAlgorithm] = None,
            moalgo: Optional[BinaryOptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedSelectionProtocol.

        Parameters
        ----------
        ntrait : Integral
            Number of traits to expect from breeding value matrix inputs.

        nself : Integral
            Number of selfing generations post-cross before progenies are derived. 
            Used in variance calculations.
        
        upper_percentile : Real
            Upper percentile of progenies to select.
        
        vmatfcty : GeneticVarianceMatrixFactory
            Variance matrix factory to calculate progeny genetic variances.
        
        gmapfn : GeneticMapFunction
            Genetic map function used in the calculation of progeny genetic variance matrices.
        
        unique_parents : bool
            Whether parents should be unique.

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

        soalgo : BinaryOptimizationAlgorithm, None
            Single-objective optimization algorithm.

            If ``soalgo`` is ``None``, then use a default single-objective 
            optimization algorithm.

        moalgo : BinaryOptimizationAlgorithm, None
            Multi-objective opimization algorithm.

            If ``moalgo`` is ``None``, then use a default multi-objective 
            optimization algorithm.

        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        # make assignments from Mixin class first
        self.ntrait = ntrait
        self.nself = nself
        self.upper_percentile = upper_percentile
        self.vmatfcty = vmatfcty
        self.gmapfn = gmapfn
        self.unique_parents = unique_parents
        # make assignments from BinaryMateSelectionProtocol second
        super(UsefulnessCriterionBinarySelection, self).__init__(
            ncross = ncross,
            nparent = nparent,
            nmating = nmating,
            nprogeny = nprogeny,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            ndset_wt = ndset_wt,
            ndset_trans = ndset_trans,
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
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
        ) -> BinaryMateSelectionProblem:
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
        out : BinaryMateSelectionProblem
            An optimization problem definition.
        """
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # get the cross map (inefficient)
        xmap = UsefulnessCriterionBinaryMateSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, len(xmap))
        decn_space_upper = numpy.repeat(1, len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # get the median number of mating from the mating property
        nmating_median = round(numpy.median(self.nmating))
        nprogeny_median = round(numpy.median(self.nprogeny))

        # construct problem
        prob = UsefulnessCriterionBinaryMateSelectionProblem.from_pgmat_gpmod_xmap(
            nparent = self.nparent, 
            ncross = nmating_median, 
            nprogeny = nprogeny_median, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            xmap = xmap,
            ndecn = len(xmap),
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = self.nobj,
            obj_wt = self.obj_wt,
            obj_trans = self.obj_trans,
            obj_trans_kwargs = self.obj_trans_kwargs,
            nineqcv = self.nineqcv,
            ineqcv_wt = self.ineqcv_wt,
            ineqcv_trans = self.ineqcv_trans,
            ineqcv_trans_kwargs = self.ineqcv_trans_kwargs,
            neqcv = self.neqcv,
            eqcv_wt = self.eqcv_wt,
            eqcv_trans = self.eqcv_trans,
            eqcv_trans_kwargs = self.eqcv_trans_kwargs
        )

        return prob

    ################ Single Objective Solve ################
    # inherit sosolve() from BinaryMateSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from BinaryMateSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from BinaryMateSelectionProtocol

class UsefulnessCriterionIntegerSelection(UsefulnessCriterionSelectionMixin,IntegerMateSelectionProtocol):
    """
    Class defining Usefulness Criterion (UC) Selection for a integer search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool,
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
            soalgo: Optional[IntegerOptimizationAlgorithm] = None,
            moalgo: Optional[IntegerOptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedSelectionProtocol.

        Parameters
        ----------
        ntrait : Integral
            Number of traits to expect from breeding value matrix inputs.

        nself : Integral
            Number of selfing generations post-cross before progenies are derived. 
            Used in variance calculations.
        
        upper_percentile : Real
            Upper percentile of progenies to select.
        
        vmatfcty : GeneticVarianceMatrixFactory
            Variance matrix factory to calculate progeny genetic variances.
        
        gmapfn : GeneticMapFunction
            Genetic map function used in the calculation of progeny genetic variance matrices.
        
        unique_parents : bool
            Whether parents should be unique.

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

        soalgo : IntegerOptimizationAlgorithm, None
            Single-objective optimization algorithm.

            If ``soalgo`` is ``None``, then use a default single-objective 
            optimization algorithm.

        moalgo : IntegerOptimizationAlgorithm, None
            Multi-objective opimization algorithm.

            If ``moalgo`` is ``None``, then use a default multi-objective 
            optimization algorithm.

        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        # make assignments from Mixin class first
        self.ntrait = ntrait
        self.nself = nself
        self.upper_percentile = upper_percentile
        self.vmatfcty = vmatfcty
        self.gmapfn = gmapfn
        self.unique_parents = unique_parents
        # make assignments from IntegerMateSelectionProtocol second
        super(UsefulnessCriterionIntegerSelection, self).__init__(
            ncross = ncross,
            nparent = nparent,
            nmating = nmating,
            nprogeny = nprogeny,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            ndset_wt = ndset_wt,
            ndset_trans = ndset_trans,
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
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
        ) -> IntegerMateSelectionProblem:
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
        out : IntegerMateSelectionProblem
            An optimization problem definition.
        """
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # get the cross map (inefficient)
        xmap = UsefulnessCriterionIntegerMateSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, len(xmap))
        decn_space_upper = numpy.repeat(self.ncross * self.nparent * self.nmating, len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # get the median number of mating from the mating property
        nmating_median = round(numpy.median(self.nmating))
        nprogeny_median = round(numpy.median(self.nprogeny))

        # construct problem
        prob = UsefulnessCriterionIntegerMateSelectionProblem.from_pgmat_gpmod_xmap(
            nparent = self.nparent, 
            ncross = nmating_median, 
            nprogeny = nprogeny_median, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            xmap = xmap,
            ndecn = len(xmap),
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = self.nobj,
            obj_wt = self.obj_wt,
            obj_trans = self.obj_trans,
            obj_trans_kwargs = self.obj_trans_kwargs,
            nineqcv = self.nineqcv,
            ineqcv_wt = self.ineqcv_wt,
            ineqcv_trans = self.ineqcv_trans,
            ineqcv_trans_kwargs = self.ineqcv_trans_kwargs,
            neqcv = self.neqcv,
            eqcv_wt = self.eqcv_wt,
            eqcv_trans = self.eqcv_trans,
            eqcv_trans_kwargs = self.eqcv_trans_kwargs
        )

        return prob

    ################ Single Objective Solve ################
    # inherit sosolve() from IntegerMateSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from IntegerMateSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from IntegerMateSelectionProtocol

class UsefulnessCriterionRealSelection(UsefulnessCriterionSelectionMixin,RealMateSelectionProtocol):
    """
    Class defining Usefulness Criterion (UC) Selection for real search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool,
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
            soalgo: Optional[RealOptimizationAlgorithm] = None,
            moalgo: Optional[RealOptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedSelectionProtocol.

        Parameters
        ----------
        ntrait : Integral
            Number of traits to expect from breeding value matrix inputs.

        nself : Integral
            Number of selfing generations post-cross before progenies are derived. 
            Used in variance calculations.
        
        upper_percentile : Real
            Upper percentile of progenies to select.
        
        vmatfcty : GeneticVarianceMatrixFactory
            Variance matrix factory to calculate progeny genetic variances.
        
        gmapfn : GeneticMapFunction
            Genetic map function used in the calculation of progeny genetic variance matrices.
        
        unique_parents : bool
            Whether parents should be unique.

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

        soalgo : RealOptimizationAlgorithm, None
            Single-objective optimization algorithm.

            If ``soalgo`` is ``None``, then use a default single-objective 
            optimization algorithm.

        moalgo : RealOptimizationAlgorithm, None
            Multi-objective opimization algorithm.

            If ``moalgo`` is ``None``, then use a default multi-objective 
            optimization algorithm.

        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        # make assignments from Mixin class first
        self.ntrait = ntrait
        self.nself = nself
        self.upper_percentile = upper_percentile
        self.vmatfcty = vmatfcty
        self.gmapfn = gmapfn
        self.unique_parents = unique_parents
        # make assignments from RealMateSelectionProtocol second
        super(UsefulnessCriterionRealSelection, self).__init__(
            ncross = ncross,
            nparent = nparent,
            nmating = nmating,
            nprogeny = nprogeny,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            ndset_wt = ndset_wt,
            ndset_trans = ndset_trans,
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
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
        ) -> RealMateSelectionProblem:
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
        out : RealMateSelectionProblem
            An optimization problem definition.
        """
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # get the cross map (inefficient)
        xmap = UsefulnessCriterionRealMateSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0.0, len(xmap))
        decn_space_upper = numpy.repeat(1.0, len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # get the median number of mating from the mating property
        nmating_median = round(numpy.median(self.nmating))
        nprogeny_median = round(numpy.median(self.nprogeny))

        # construct problem
        prob = UsefulnessCriterionRealMateSelectionProblem.from_pgmat_gpmod_xmap(
            nparent = self.nparent, 
            ncross = nmating_median, 
            nprogeny = nprogeny_median, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            xmap = xmap,
            ndecn = len(xmap),
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = self.nobj,
            obj_wt = self.obj_wt,
            obj_trans = self.obj_trans,
            obj_trans_kwargs = self.obj_trans_kwargs,
            nineqcv = self.nineqcv,
            ineqcv_wt = self.ineqcv_wt,
            ineqcv_trans = self.ineqcv_trans,
            ineqcv_trans_kwargs = self.ineqcv_trans_kwargs,
            neqcv = self.neqcv,
            eqcv_wt = self.eqcv_wt,
            eqcv_trans = self.eqcv_trans,
            eqcv_trans_kwargs = self.eqcv_trans_kwargs
        )

        return prob

    ################ Single Objective Solve ################
    # inherit sosolve() from RealMateSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from RealMateSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from RealMateSelectionProtocol

class UsefulnessCriterionSubsetSelection(UsefulnessCriterionSelectionMixin,SubsetMateSelectionProtocol):
    """
    Class defining Usefulness Criterion (UC) Selection for subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ntrait: Integral,
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool,
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
            soalgo: Optional[SubsetOptimizationAlgorithm] = None,
            moalgo: Optional[SubsetOptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedSelectionProtocol.

        Parameters
        ----------
        ntrait : Integral
            Number of traits to expect from breeding value matrix inputs.

        nself : Integral
            Number of selfing generations post-cross before progenies are derived. 
            Used in variance calculations.
        
        upper_percentile : Real
            Upper percentile of progenies to select.
        
        vmatfcty : GeneticVarianceMatrixFactory
            Variance matrix factory to calculate progeny genetic variances.
        
        gmapfn : GeneticMapFunction
            Genetic map function used in the calculation of progeny genetic variance matrices.
        
        unique_parents : bool
            Whether parents should be unique.

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

        soalgo : SubsetOptimizationAlgorithm, None
            Single-objective optimization algorithm.

            If ``soalgo`` is ``None``, then use a default single-objective 
            optimization algorithm.

        moalgo : SubsetOptimizationAlgorithm, None
            Multi-objective opimization algorithm.

            If ``moalgo`` is ``None``, then use a default multi-objective 
            optimization algorithm.

        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        # make assignments from Mixin class first
        self.ntrait = ntrait
        self.nself = nself
        self.upper_percentile = upper_percentile
        self.vmatfcty = vmatfcty
        self.gmapfn = gmapfn
        self.unique_parents = unique_parents
        # make assignments from SubsetMateSelectionProtocol second
        super(UsefulnessCriterionSubsetSelection, self).__init__(
            ncross = ncross,
            nparent = nparent,
            nmating = nmating,
            nprogeny = nprogeny,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            ndset_wt = ndset_wt,
            ndset_trans = ndset_trans,
            ndset_trans_kwargs = ndset_trans_kwargs,
            rng = rng,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
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
        ) -> SubsetMateSelectionProblem:
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
        out : SubsetMateSelectionProblem
            An optimization problem definition.
        """
        # type checks
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_GenomicModel(gpmod, "gpmod")
        
        # get the cross map
        xmap = UsefulnessCriterionSubsetMateSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space = numpy.arange(len(xmap))
        decn_space_lower = numpy.repeat(0, self.ncross)
        decn_space_upper = numpy.repeat(len(xmap)-1, self.ncross)

        # get the median number of mating from the mating property
        nmating_median = round(numpy.median(self.nmating))
        nprogeny_median = round(numpy.median(self.nprogeny))

        # construct problem
        prob = UsefulnessCriterionSubsetMateSelectionProblem.from_pgmat_gpmod_xmap(
            nparent = self.nparent, 
            ncross = nmating_median,
            nprogeny = nprogeny_median, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            xmap = xmap,
            ndecn = self.ncross,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = self.nobj,
            obj_wt = self.obj_wt,
            obj_trans = self.obj_trans,
            obj_trans_kwargs = self.obj_trans_kwargs,
            nineqcv = self.nineqcv,
            ineqcv_wt = self.ineqcv_wt,
            ineqcv_trans = self.ineqcv_trans,
            ineqcv_trans_kwargs = self.ineqcv_trans_kwargs,
            neqcv = self.neqcv,
            eqcv_wt = self.eqcv_wt,
            eqcv_trans = self.eqcv_trans,
            eqcv_trans_kwargs = self.eqcv_trans_kwargs
        )

        return prob

    ################ Single Objective Solve ################
    # inherit sosolve() from SubsetMateSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from SubsetMateSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from SubsetMateSelectionProtocol
