"""
Module defining mean expected heterozygosity selection protocols.
"""

# list of all public objects in this module
__all__ = [
    "L2NormGenomicSelectionMixin",
    "L2NormGenomicBinarySelection",
    "L2NormGenomicIntegerSelection",
    "L2NormGenomicRealSelection",
    "L2NormGenomicSubsetSelection",
]

# imports
from abc import ABCMeta
from numbers import Integral
from numbers import Real
from typing import Callable
from typing import Optional
from typing import Union

import numpy
from numpy.random import Generator
from numpy.random import RandomState
import pandas

from pybrops.breed.prot.sel.BinarySelectionProtocol import BinarySelectionProtocol
from pybrops.breed.prot.sel.IntegerSelectionProtocol import IntegerSelectionProtocol
from pybrops.breed.prot.sel.RealSelectionProtocol import RealSelectionProtocol
from pybrops.breed.prot.sel.SubsetSelectionProtocol import SubsetSelectionProtocol
from pybrops.breed.prot.sel.prob.BinarySelectionProblem import BinarySelectionProblem
from pybrops.breed.prot.sel.prob.IntegerSelectionProblem import IntegerSelectionProblem
from pybrops.breed.prot.sel.prob.RealSelectionProblem import RealSelectionProblem
from pybrops.breed.prot.sel.prob.L2NormGenomicSelectionProblem import L2NormGenomicBinarySelectionProblem
from pybrops.breed.prot.sel.prob.L2NormGenomicSelectionProblem import L2NormGenomicIntegerSelectionProblem
from pybrops.breed.prot.sel.prob.L2NormGenomicSelectionProblem import L2NormGenomicRealSelectionProblem
from pybrops.breed.prot.sel.prob.L2NormGenomicSelectionProblem import L2NormGenomicSubsetSelectionProblem
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.opt.algo.BinaryOptimizationAlgorithm import BinaryOptimizationAlgorithm
from pybrops.opt.algo.IntegerOptimizationAlgorithm import IntegerOptimizationAlgorithm
from pybrops.opt.algo.RealOptimizationAlgorithm import RealOptimizationAlgorithm
from pybrops.opt.algo.SubsetOptimizationAlgorithm import SubsetOptimizationAlgorithm
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import CoancestryMatrixFactory
from pybrops.popgen.cmat.fcty.CoancestryMatrixFactory import check_is_CoancestryMatrixFactory
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class L2NormGenomicSelectionMixin(
        metaclass = ABCMeta,
    ):
    """
    Semi-abstract class for L2-Norm Genomic Selection (L2GS) with constraints.
    """

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################
    @property
    def cmatfcty(self) -> CoancestryMatrixFactory:
        """Coancestry matrix factory used to create coancestry matrices."""
        return self._cmatfcty
    @cmatfcty.setter
    def cmatfcty(self, value: CoancestryMatrixFactory) -> None:
        """Set coancestry matrix factory."""
        check_is_CoancestryMatrixFactory(value, "cmatfcty")
        self._cmatfcty = value

class L2NormGenomicBinarySelection(L2NormGenomicSelectionMixin,BinarySelectionProtocol):
    """
    L2-Norm Genomic Selection (L2GS) in a subset search space.
    """
    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            cmatfcty: CoancestryMatrixFactory,
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
        Constructor for the concrete class L2NormGenomicBinarySelection.

        Parameters
        ----------
        cmatfcty : CoancestryMatrixFactory
            Coancestry matrix factory used to create coancestry matrices.
        
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
        self.cmatfcty = cmatfcty
        # make assignments from BinarySelectionProtocol second
        super(L2NormGenomicBinarySelection, self).__init__(
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
        ) -> BinarySelectionProblem:
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
        out : BinarySelectionProblem
            An optimization problem definition.
        """
        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, ntaxa)
        decn_space_upper = numpy.repeat(1, ntaxa)
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = L2NormGenomicBinarySelectionProblem.from_gmat(
            gmat = gmat,
            cmatfcty = self.cmatfcty,
            ndecn = self.nparent,
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
    # inherit sosolve() from BinarySelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from BinarySelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from BinarySelectionProtocol

class L2NormGenomicIntegerSelection(L2NormGenomicSelectionMixin,IntegerSelectionProtocol):
    """
    L2-Norm Genomic Selection (L2GS) in an integer search space.
    """
    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            cmatfcty: CoancestryMatrixFactory,
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
        Constructor for the concrete class L2NormGenomicIntegerSelection.

        Parameters
        ----------
        cmatfcty : CoancestryMatrixFactory
            Coancestry matrix factory used to create coancestry matrices.
        
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
        self.cmatfcty = cmatfcty
        # make assignments from IntegerSelectionProtocol second
        super(L2NormGenomicIntegerSelection, self).__init__(
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
        ) -> IntegerSelectionProblem:
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
        out : IntegerSelectionProblem
            An optimization problem definition.
        """
        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, ntaxa)
        decn_space_upper = numpy.repeat(ntaxa, ntaxa)
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = L2NormGenomicIntegerSelectionProblem.from_gmat(
            gmat = gmat,
            cmatfcty = self.cmatfcty,
            ndecn = self.nparent,
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
    # inherit sosolve() from IntegerSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from IntegerSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from IntegerSelectionProtocol

class L2NormGenomicRealSelection(L2NormGenomicSelectionMixin,RealSelectionProtocol):
    """
    L2-Norm Genomic Selection (L2GS) in a real search space.
    """
    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            cmatfcty: CoancestryMatrixFactory,
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
        Constructor for the concrete class L2NormGenomicRealSelection.

        Parameters
        ----------
        cmatfcty : CoancestryMatrixFactory
            Coancestry matrix factory used to create coancestry matrices.
        
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
        self.cmatfcty = cmatfcty
        # make assignments from RealSelectionProtocol second
        super(L2NormGenomicRealSelection, self).__init__(
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
        ) -> RealSelectionProblem:
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
        out : RealSelectionProblem
            An optimization problem definition.
        """
        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space_lower = numpy.repeat(0.0, ntaxa)
        decn_space_upper = numpy.repeat(1.0, ntaxa)
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = L2NormGenomicRealSelectionProblem.from_gmat(
            gmat = gmat,
            cmatfcty = self.cmatfcty,
            ndecn = self.nparent,
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
    # inherit sosolve() from RealSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from RealSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from RealSelectionProtocol

class L2NormGenomicSubsetSelection(L2NormGenomicSelectionMixin,SubsetSelectionProtocol):
    """
    L2-Norm Genomic Selection (L2GS) in a subset search space.
    """
    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            cmatfcty: CoancestryMatrixFactory,
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
        Constructor for the concrete class L2NormGenomicSubsetSelection.

        Parameters
        ----------
        cmatfcty : CoancestryMatrixFactory
            Coancestry matrix factory used to create coancestry matrices.
        
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
        self.cmatfcty = cmatfcty
        # make assignments from SubsetSelectionProtocol second
        super(L2NormGenomicSubsetSelection, self).__init__(
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
        ) -> SubsetSelectionProblem:
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
        out : SubsetSelectionProblem
            An optimization problem definition.
        """
        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space = numpy.arange(ntaxa)
        decn_space_lower = numpy.repeat(0, self.nparent)
        decn_space_upper = numpy.repeat(ntaxa-1, self.nparent)

        # construct problem
        prob = L2NormGenomicSubsetSelectionProblem.from_gmat(
            gmat = gmat,
            cmatfcty = self.cmatfcty,
            ndecn = self.nparent,
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
    # inherit sosolve() from SubsetSelectionProtocol

    ################ Multi Objective Solve #################
    # inherit mosolve() from SubsetSelectionProtocol

    ################# Selection Functions ##################
    # inherit select() from SubsetSelectionProtocol
