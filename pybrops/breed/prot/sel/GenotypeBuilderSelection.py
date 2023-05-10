"""
Module defining protocols for Genotype Builder (GB) selection.
"""

__all__ = [

]

from abc import ABCMeta

from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol


class GenotypeBuilderBaseSelection(SelectionProtocol,metaclass=ABCMeta):
    """
    docstring for GenotypeBuilderBaseSelection.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for GenotypeBuilderBaseSelection.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(GenotypeBuilderBaseSelection, self).__init__(**kwargs)

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################
