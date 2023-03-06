"""
Module defining interfaces and error checking routines for selection protocols.
"""

from numbers import Integral, Real
from typing import Any, Callable, Optional, Union
import numpy
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

class ConstrainedSelectionProtocol:
    """
    Abstract class defining interfaces for selection protocols.

    The purpose of this abstract class is to define functionality for:
        1) Selection (both single- and multi-objective) of genotypes.
        2) Construction of fast objective functions.
        3) Mapping of the Pareto frontier for the selection protocol.
        4) Access to static objective functions for the selection protocol.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedSelectionProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(ConstrainedSelectionProtocol, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def nobj(self) -> Integral:
        """Number of objectives. Equivalent to length of vector output by ``objfn_trans``."""
        raise NotImplementedError("property is abstract")
    @nobj.setter
    def nobj(self, value: Integral) -> None:
        """Set number of objectives."""
        raise NotImplementedError("property is abstract")
    @nobj.deleter
    def nobj(self) -> None:
        """Delete number of objectives."""
        raise NotImplementedError("property is abstract")
    
    @property
    def objfn_trans(self) -> Callable[[numpy.ndarray,dict],numpy.ndarray]:
        """Function for transforming latent vector into objective function value(s)."""
        raise NotImplementedError("property is abstract")
    @objfn_trans.setter
    def objfn_trans(self, value: Callable[[numpy.ndarray,dict],numpy.ndarray]) -> None:
        """Set function for transforming latent vector into objective function value(s)."""
        raise NotImplementedError("property is abstract")
    @objfn_trans.deleter
    def objfn_trans(self) -> None:
        """Delete function for transforming latent vector into objective function value(s)."""
        raise NotImplementedError("property is abstract")
    
    @property
    def objfn_trans_kwargs(self) -> dict:
        """Function keyword arguments for transforming latent vector into objective function value(s)."""
        raise NotImplementedError("property is abstract")
    @objfn_trans_kwargs.setter
    def objfn_trans_kwargs(self, value: dict) -> None:
        """Set function keyword arguments for transforming latent vector into objective function value(s)."""
        raise NotImplementedError("property is abstract")
    @objfn_trans_kwargs.deleter
    def objfn_trans_kwargs(self) -> None:
        """Delete function keyword arguments for transforming latent vector into objective function value(s)."""
        raise NotImplementedError("property is abstract")

    @property
    def nineqcv(self) -> Integral:
        """Number of inequality constraint violations. Equivalent to length of vector output by ``ineqcvfn_trans``."""
        raise NotImplementedError("property is abstract")
    @nineqcv.setter
    def nineqcv(self, value: Integral) -> None:
        """Set number of inequality constraint violations."""
        raise NotImplementedError("property is abstract")
    @nineqcv.deleter
    def nineqcv(self) -> None:
        """Delete number of inequality constraint violations."""
        raise NotImplementedError("property is abstract")

    @property
    def ineqcvfn_trans(self) -> Callable[[numpy.ndarray,dict],numpy.ndarray]:
        """Function for transforming latent vector into inequality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")
    @ineqcvfn_trans.setter
    def ineqcvfn_trans(self, value: Callable[[numpy.ndarray,dict],numpy.ndarray]) -> None:
        """Set function for transforming latent vector into inequality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")
    @ineqcvfn_trans.deleter
    def ineqcvfn_trans(self) -> None:
        """Delete function for transforming latent vector into inequality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")
    
    @property
    def ineqcvfn_trans_kwargs(self) -> dict:
        """Function keyword arguments for transforming latent vector into inequality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")
    @ineqcvfn_trans_kwargs.setter
    def ineqcvfn_trans_kwargs(self, value: dict) -> None:
        """Set function keyword arguments for transforming latent vector into inequality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")
    @ineqcvfn_trans_kwargs.deleter
    def ineqcvfn_trans_kwargs(self) -> None:
        """Delete function keyword arguments for transforming latent vector into inequality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")

    @property
    def neqcv(self) -> Integral:
        """Number of equality constraint violations. Equivalent to length of vector output by ``eqcvfn_trans``."""
        raise NotImplementedError("property is abstract")
    @neqcv.setter
    def neqcv(self, value: Integral) -> None:
        """Set number of equality constraint violations."""
        raise NotImplementedError("property is abstract")
    @neqcv.deleter
    def neqcv(self) -> None:
        """Delete number of equality constraint violations."""
        raise NotImplementedError("property is abstract")
    
    @property
    def eqcvfn_trans(self) -> Callable[[numpy.ndarray,dict],numpy.ndarray]:
        """Function for transforming latent vector into equality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")
    @eqcvfn_trans.setter
    def eqcvfn_trans(self, value: Callable[[numpy.ndarray,dict],numpy.ndarray]) -> None:
        """Set function for transforming latent vector into equality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")
    @eqcvfn_trans.deleter
    def eqcvfn_trans(self) -> None:
        """Delete function for transforming latent vector into equality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")
    
    @property
    def eqcvfn_trans_kwargs(self) -> dict:
        """Function keyword arguments for transforming latent vector into equality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")
    @eqcvfn_trans_kwargs.setter
    def eqcvfn_trans_kwargs(self, value: dict) -> None:
        """Set function keyword arguments for transforming latent vector into equality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")
    @eqcvfn_trans_kwargs.deleter
    def eqcvfn_trans_kwargs(self) -> None:
        """Delete function keyword arguments for transforming latent vector into equality constraint violation function value(s)."""
        raise NotImplementedError("property is abstract")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Function evaluation and construction #########
    def encodefn(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            **kwargs: dict
        ) -> Callable:
        """
        Return a latent space encoding function constructed from inputs for fast
        calling. The purpose of this function is to create function that converts
        a candidate solution vector into a latent space vector. This latent space
        vector contains all information needed to construct objective functions,
        inequality constraints, and equality constraints.

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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Callable
            A selection objective function for the specified problem.
        """
        raise NotImplementedError("method is abstract")

    def evalfn(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            **kwargs: dict
        ) -> Callable:
        """
        Return an evaluation function constructed from inputs for fast calling. 
        The purpose of this function is to create function that converts a
        candidate solution vector into a tuple of scalars and/or vectors. The 
        returned tuple contains three elements: an objective function scalar or
        vector, an inequality constraint violation scalar or vector, and an 
        equality constraint violation scalar or vector.

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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Callable
            A selection objective function for the specified problem.
        """
        raise NotImplementedError("method is abstract")

    def objfn(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            **kwargs: dict
        ) -> Callable:
        """
        Return an objective function constructed from inputs for fast calling.
        The returned function converts a candidate solution vector into a scalar
        or vector of objective value(s). This returned constructed function does
        not return any inequality or equality constraint violation value(s).
        To obtain multiple value types, use the ``evalfn`` function.

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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Callable
            A candidate solution objective function for the specified problem.
        """
        raise NotImplementedError("method is abstract")

    def ineqcvfn(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            **kwargs: dict
        ) -> Callable:
        """
        Return an inequality constraint violation function constructed from 
        inputs for fast calling. The returned function converts a candidate 
        solution vector into a scalar or vector of inequality constraint 
        violation value(s). This returned constructed function does not return
        any objective or equality constraint violation value(s).
        To obtain multiple value types, use the ``evalfn`` function.

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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Callable
            A candidate solution inequality constraint violation function for 
            the specified problem.
        """
        raise NotImplementedError("method is abstract")

    def eqcvfn(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            **kwargs: dict
        ) -> Callable:
        """
        Return an equality constraint violation function constructed from inputs 
        for fast calling. The returned function converts a candidate solution 
        vector into a scalar or vector of equality constraint violation value(s). 
        This returned constructed function does not return any objective or 
        equality constraint violation value(s).
        To obtain multiple value types, use the ``evalfn`` function.

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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Callable
            A candidate solution equality constraint violation function for the 
            specified problem.
        """
        raise NotImplementedError("method is abstract")

    ################# Selection Functions ##################
    def select(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict], 
            **kwargs: dict
        ) -> tuple:
        """
        Select individuals for breeding.

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
        miscout : dict, None
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
        raise NotImplementedError("method is abstract")

    ############## Pareto Frontier Functions ###############
    def pareto(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict], 
            **kwargs: dict
        ) -> tuple:
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
        miscout : dict, None
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
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def encodefn_static(sel: numpy.ndarray, *args: tuple, **kwargs: dict) -> numpy.ndarray:
        """
        Encoding function for the selection protocol.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection vector of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals.
        args : tuple
            Additional positional arguments.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of encoded space values of shape ``(o,)`` or a scalar.

            Where:

            - ``o`` is the number of objectives.
        """
        raise NotImplementedError("static method is abstract")
    
    @staticmethod
    def evalfn_static(sel: numpy.ndarray, *args: tuple, **kwargs: dict) -> tuple:
        """
        Evaluation function for the selection protocol.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection vector of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals.
        args : tuple
            Additional positional arguments.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple ``(obj, ineqcv, eqcv)`` that contains objective function
            evaluations (``obj``), inequality constraint violation values
            (``ineqcv``), and equality constraint violation values (``eqcv``).
        """
        raise NotImplementedError("static method is abstract")

    @staticmethod
    def objfn_static(sel: numpy.ndarray, *args: tuple, **kwargs: dict) -> Union[Real,numpy.ndarray]:
        """
        Objective function for the selection protocol.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection vector of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals.
        args : tuple
            Additional positional arguments.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray, Real
            An array of objective function values of shape ``(o,)`` or a scalar.

            Where:

            - ``o`` is the number of objectives.
        """
        raise NotImplementedError("static method is abstract")

    @staticmethod
    def ineqcv_static(sel: numpy.ndarray, *args: tuple, **kwargs: dict) -> Union[Real,numpy.ndarray]:
        """
        Inequality constraint violation function for the selection protocol.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection vector of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals.
        args : tuple
            Additional positional arguments.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray, Real
            An array of objective function values of shape ``(o,)`` or a scalar.

            Where:

            - ``o`` is the number of objectives.
        """
        raise NotImplementedError("static method is abstract")

    @staticmethod
    def eqcv_static(sel: numpy.ndarray, *args: tuple, **kwargs: dict) -> Union[Real,numpy.ndarray]:
        """
        Equality constraint violation function for the selection protocol.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection vector of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals.
        args : tuple
            Additional positional arguments.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray, Real
            An array of objective function values of shape ``(o,)`` or a scalar.

            Where:

            - ``o`` is the number of objectives.
        """
        raise NotImplementedError("static method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_NewSelectionProtocol(v: Any) -> bool:
    return isinstance(v, ConstrainedSelectionProtocol)

def check_is_NewSelectionProtocol(v: Any, vname: str) -> None:
    if not isinstance(v, ConstrainedSelectionProtocol):
        raise TypeError("variable '{0}' must be a ConstrainedSelectionProtocol".format(vname))
