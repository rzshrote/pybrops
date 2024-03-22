"""
Module defining a general class for selection protocols.
"""

# list of all public objects in this module
__all__ = [
    "MateSelectionProtocol",
    "check_is_MateSelectionProtocol",
]

# imports
from abc import ABCMeta
from abc import abstractmethod
from numbers import Integral
from typing import Optional

import pandas

from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.breed.prot.sel.cfg.MateSelectionConfiguration import MateSelectionConfiguration
from pybrops.breed.prot.sel.prob.MateSelectionProblem import MateSelectionProblem
from pybrops.breed.prot.sel.soln.MateSelectionSolution import MateSelectionSolution
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class MateSelectionProtocol(
        SelectionProtocol,
        metaclass = ABCMeta,
    ):
    """
    A semi-abstract class implementing several key properties common to most, 
    if not all, constrained selection protocols.
    """

    ########################## Special Object Methods ##########################
    # inherit from SelectionProtocol

    ############################ Object Properties #############################
    # inherit from SelectionProtocol

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
        ) -> MateSelectionProblem:
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
        ) -> MateSelectionSolution:
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
        ) -> MateSelectionSolution:
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
        ) -> MateSelectionConfiguration:
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
def check_is_MateSelectionProtocol(v: object, vname: str) -> None:
    """
    Check if object is of type ConstrainedMateSelectionProtocol, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, MateSelectionProtocol):
        raise TypeError("'{0}' must be of type ConstrainedMateSelectionProtocol.".format(vname))
