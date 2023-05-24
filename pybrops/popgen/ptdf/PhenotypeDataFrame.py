"""
Module defining phenotype dataframe interfaces and associated error checking routines.
"""

from typing import Any
from pybrops.core.df.DataFrame import DataFrame

class PhenotypeDataFrame(DataFrame):
    """
    An abstract class for phenotype data frame objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Data frame column analysis and effect types
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class PhenotypeDataFrame.
        """
        super(PhenotypeDataFrame, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    # TODO: maybe eliminate these. it seems like this information should go in other modules
    @property
    def col_analysis_type(self) -> Any:
        """Analysis variable type array."""
        raise NotImplementedError("property is abstract")
    @col_analysis_type.setter
    def col_analysis_type(self, value: Any) -> None:
        """Set analysis variable type array"""
        raise NotImplementedError("property is abstract")
    @col_analysis_type.deleter
    def col_analysis_type(self) -> None:
        """Delete analysis variable type array"""
        raise NotImplementedError("property is abstract")

    # TODO: maybe eliminate these. it seems like this information should go in other modules
    @property
    def col_analysis_effect(self) -> Any:
        """Analysis variable effect type {'response','fixed','random',None} array."""
        raise NotImplementedError("property is abstract")
    @col_analysis_effect.setter
    def col_analysis_effect(self, value: Any) -> None:
        """Set analysis variable effect type array"""
        raise NotImplementedError("property is abstract")
    @col_analysis_effect.deleter
    def col_analysis_effect(self) -> None:
        """Delete analysis variable effect type array"""
        raise NotImplementedError("property is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhenotypeDataFrame(v: object) -> bool:
    """
    Determine whether an object is a PhenotypeDataFrame.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a PhenotypeDataFrame object instance.
    """
    return isinstance(v, PhenotypeDataFrame)

def check_is_PhenotypeDataFrame(v: object, vname: str) -> None:
    """
    Check if object is of type PhenotypeDataFrame. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, PhenotypeDataFrame):
        raise TypeError("variable '{0}' must be a PhenotypeDataFrame".format(vname))
