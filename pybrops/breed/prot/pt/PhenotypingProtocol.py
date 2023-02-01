"""
Module defining interfaces and associated protocols for phenotyping protocols.
"""

from typing import Any


class PhenotypingProtocol:
    """
    Abstract class defining interfaces for phenotyping protocols.

    The purpose of this abstract class is to provide functionality for:
        1) Genomic model metadata.
        2) Phenotype simulation.
        3) Manipulation and setting of environmental variance metadata.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict) -> None:
        """
        Constructor for the abstract class PhenotypingProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(PhenotypingProtocol, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############### Genomic Model Properties ###############
    def gpmod():
        doc = "Genomic prediction model."
        def fget(self):
            """Get genomic prediction model"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set genomic prediction model"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete genomic prediction model"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    gpmod = property(**gpmod())

    ################ Stochastic Parameters #################
    def var_err():
        doc = "Error variance for each trait."
        def fget(self):
            """Get error variance"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set error variance"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete error variance"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    var_err = property(**var_err())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def phenotype(self, pgmat, miscout, **kwargs: dict):
        """
        Phenotype a set of genotypes using a genomic prediction model.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes of the individuals to phenotype.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhenotypeDataFrame
            A PhenotypeDataFrame containing phenotypes for individuals.
        """
        raise NotImplementedError("method is abstract")

    def set_h2(self, h2, pgmat, **kwargs: dict):
        """
        Set the narrow sense heritability for environments.

        Parameters
        ----------
        h2 : float, numpy.ndarray
            Narrow sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise NotImplementedError("method is abstract")

    def set_H2(self, H2, pgmat, **kwargs: dict):
        """
        Set the broad sense heritability for environments.

        Parameters
        ----------
        H2 : float, numpy.ndarray
            Broad sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhenotypingProtocol(v: Any) -> bool:
    """
    Determine whether an object is a PhenotypingProtocol.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a PhenotypingProtocol object instance.
    """
    return isinstance(v, PhenotypingProtocol)

def check_is_PhenotypingProtocol(v: Any, varname: str) -> None:
    """
    Check if object is of type PhenotypingProtocol. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, PhenotypingProtocol):
        raise TypeError("'%s' must be a PhenotypingProtocol." % varname)
