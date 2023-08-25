"""
Module defining interfaces and error checking routines for selection protocols.
"""

from typing import Optional
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class UnconstrainedSelectionProtocol:
    """
    Abstract class defining interfaces for selection protocols.

    The purpose of this abstract class is to define functionality for:
        1) Selection (both single- and multi-objective) of genotypes.
        2) Construction of fast objective functions.
        3) Mapping of the Pareto frontier for the selection protocol.
        4) Access to static objective functions for the selection protocol.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class SelectionProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(UnconstrainedSelectionProtocol, self).__init__()

    ############################ Object Properties #############################

    ############################## Object Methods ##############################
    def select(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: int, 
            t_max: int, 
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

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return an objective function constructed from inputs for fast calling.

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
        out : function
            A selection objective function for the specified problem.
        """
        raise NotImplementedError("method is abstract")

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return a vectorized objective function constructed from inputs for fast
        calling.

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
        out : function
            A vectorized selection objective function for the specified problem.
        """
        raise NotImplementedError("method is abstract")

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout, **kwargs: dict):
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
    def objfn_static(sel, *args):
        """
        Objective function for the selection protocol.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection vector of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals.
        args : tuple
            Additional arguments.

        Returns
        -------
        out : numpy.ndarray, scalar
            An array of objective function values of shape ``(o,)`` or a scalar.

            Where:

            - ``o`` is the number of objectives.
        """
        raise NotImplementedError("method is abstract")

    @staticmethod
    def objfn_vec_static(sel, *args):
        """
        Vectorized objective function for the selection protocol.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection vector of shape ``(j,k)``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``k`` is the number of individuals.
        args : tuple
            Additional arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of objective function values of shape ``(j,o)`` or
            ``(j,)``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``o`` is the number of objectives.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_SelectionProtocol(v: object, vname: str) -> None:
    if not isinstance(v, UnconstrainedSelectionProtocol):
        raise TypeError("variable '{0}' must be a SelectionProtocol".format(vname))
