"""
Module defining interfaces and associated error checking routines for matrices
with axes that are square and with taxa metadata.
"""

__all__ = [
    "SquareTaxaMatrix",
    "check_is_SquareTaxaMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from pybrops.core.mat.SquareMatrix import SquareMatrix
from pybrops.core.mat.TaxaMatrix import TaxaMatrix

class SquareTaxaMatrix(
        SquareMatrix,
        TaxaMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for matrix wrapper objects with square taxa metadata.

    The purpose of this abstract class is to merge the following interfaces:

        1. SquareMatrix
        2. TaxaMatrix
    """

    ########################## Special Object Methods ##########################

    ############## Forward numeric operators ###############
    ### __add__         inherited from ``SquareMatrix``
    ### __sub__         inherited from ``SquareMatrix``
    ### __mul__         inherited from ``SquareMatrix``
    ### __matmul__      inherited from ``SquareMatrix``
    ### __truediv__     inherited from ``SquareMatrix``
    ### __floordiv__    inherited from ``SquareMatrix``
    ### __mod__         inherited from ``SquareMatrix``
    ### __divmod__      inherited from ``SquareMatrix``
    ### __pow__         inherited from ``SquareMatrix``
    ### __lshift__      inherited from ``SquareMatrix``
    ### __rshift__      inherited from ``SquareMatrix``
    ### __and__         inherited from ``SquareMatrix``
    ### __xor__         inherited from ``SquareMatrix``
    ### __or__          inherited from ``SquareMatrix``

    ############# Backwards numeric operators ##############
    ### __radd__        inherited from ``SquareMatrix``
    ### __rsub__        inherited from ``SquareMatrix``
    ### __rmul__        inherited from ``SquareMatrix``
    ### __rmatmul__     inherited from ``SquareMatrix``
    ### __rtruediv__    inherited from ``SquareMatrix``
    ### __rfloordiv__   inherited from ``SquareMatrix``
    ### __rmod__        inherited from ``SquareMatrix``
    ### __rdivmod__     inherited from ``SquareMatrix``
    ### __rlshift__     inherited from ``SquareMatrix``
    ### __rrshift__     inherited from ``SquareMatrix``
    ### __rand__        inherited from ``SquareMatrix``
    ### __rxor__        inherited from ``SquareMatrix``
    ### __ror__         inherited from ``SquareMatrix``

    ############# Augmented numeric operators ##############
    ### __iadd__        inherited from ``SquareMatrix``
    ### __isub__        inherited from ``SquareMatrix``
    ### __imul__        inherited from ``SquareMatrix``
    ### __imatmul__     inherited from ``SquareMatrix``
    ### __itruediv__    inherited from ``SquareMatrix``
    ### __ifloordiv__   inherited from ``SquareMatrix``
    ### __imod__        inherited from ``SquareMatrix``
    ### __ipow__        inherited from ``SquareMatrix``
    ### __ilshift__     inherited from ``SquareMatrix``
    ### __irshift__     inherited from ``SquareMatrix``
    ### __iand__        inherited from ``SquareMatrix``
    ### __ixor__        inherited from ``SquareMatrix``
    ### __ior__         inherited from ``SquareMatrix``

    ################## Logical operators ###################
    ### __lt__          inherited from ``SquareMatrix``
    ### __le__          inherited from ``SquareMatrix``
    ### __eq__          inherited from ``SquareMatrix``
    ### __ne__          inherited from ``SquareMatrix``
    ### __gt__          inherited from ``SquareMatrix``
    ### __ge__          inherited from ``SquareMatrix``

    ################# Container operators ##################
    ### __len__         inherited from ``SquareMatrix``
    ### __getitem__     inherited from ``SquareMatrix``
    ### __setitem__     inherited from ``SquareMatrix``
    ### __delitem__     inherited from ``SquareMatrix``
    ### __iter__        inherited from ``SquareMatrix``

    #################### Matrix copying ####################
    ### __copy__        inherited from ``SquareMatrix``
    ### __deepcopy__    inherited from ``SquareMatrix``

    ########### Miscellaneous special functions ############
    ### __repr__        inherited from ``SquareMatrix``

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    ### mat             inherited from ``SquareMatrix``

    ############## Matrix Metadata Properties ##############
    ### mat_ndim        inherited from ``SquareMatrix``
    ### mat_shape       inherited from ``SquareMatrix``

    ############## Square Metadata Properties ##############
    ### nsquare         inherited from ``SquareMatrix``
    ### square_axes     inherited from ``SquareMatrix``
    ### square_axes_len inherited from ``SquareMatrix``

    ########### Square Taxa Metadata Properties ############
    @property
    @abstractmethod
    def nsquare_taxa(self) -> int:
        """Number of taxa axes that are square."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def square_taxa_axes(self) -> tuple:
        """Axis indices for taxa axes that are square."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def square_taxa_axes_len(self) -> tuple:
        """Axis lengths for axes that are square."""
        raise NotImplementedError("property is abstract")

    ################# Taxa Data Properites #################
    ### taxa            inherited from ``TaxaMatrix``
    ### taxa_grp        inherited from ``TaxaMatrix``

    ############### Taxa Metadata Properites ###############
    ### ntaxa           inherited from ``TaxaMatrix``
    ### taxa_axis       inherited from ``TaxaMatrix``
    ### taxa_grp_name   inherited from ``TaxaMatrix``
    ### taxa_grp_stix   inherited from ``TaxaMatrix``
    ### taxa_grp_spix   inherited from ``TaxaMatrix``
    ### taxa_grp_len    inherited from ``TaxaMatrix``

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    ### copy            inherited from ``SquareMatrix``
    ### deepcopy        inherited from ``SquareMatrix``

    ######### Matrix element copy-on-manipulation ##########
    ### adjoin          inherited from ``SquareMatrix``
    ### delete          inherited from ``SquareMatrix``
    ### insert          inherited from ``SquareMatrix``
    ### select          inherited from ``SquareMatrix``
    ### adjoin_taxa     inherited from ``TaxaMatrix``
    ### delete_taxa     inherited from ``TaxaMatrix``
    ### insert_taxa     inherited from ``TaxaMatrix``
    ### select_taxa     inherited from ``TaxaMatrix``

    ######### Matrix element in-place-manipulation #########
    ### append          inherited from ``TaxaMatrix``
    ### remove          inherited from ``TaxaMatrix``
    ### incorp          inherited from ``TaxaMatrix``
    ### append_taxa     inherited from ``TaxaMatrix``
    ### remove_taxa     inherited from ``TaxaMatrix``
    ### incorp_taxa     inherited from ``TaxaMatrix``

    ################### Sorting Methods ####################
    ### lexsort         inherited from ``TaxaMatrix``
    ### reorder         inherited from ``TaxaMatrix``
    ### sort            inherited from ``TaxaMatrix``
    ### lexsort_taxa    inherited from ``TaxaMatrix``
    ### reorder_taxa    inherited from ``TaxaMatrix``
    ### sort_taxa       inherited from ``TaxaMatrix``

    ################### Grouping Methods ###################
    ### group           inherited from ``TaxaMatrix``
    ### ungroup         inherited from ``TaxaMatrix``
    ### is_grouped      inherited from ``TaxaMatrix``
    ### group_taxa      inherited from ``TaxaMatrix``
    ### ungroup_taxa    inherited from ``TaxaMatrix``
    ### is_grouped_taxa inherited from ``TaxaMatrix``

    #################### Square Methods ####################
    ### is_square       inherited from ``SquareMatrix``

    @abstractmethod
    def is_square_taxa(
            self
        ) -> bool:
        """
        Determine whether the taxa axes lengths for the square axes are identical.

        Returns
        -------
        out : bool
            ``True`` if all square taxa axes are the same length.
            ``False`` if not all square taxa axes are the same length.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat          inherited from ``SquareMatrix``
    ### concat_taxa     inherited from ``TaxaMatrix``




################################## Utilities ###################################
def check_is_SquareTaxaMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``SquareTaxaMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, SquareTaxaMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,SquareTaxaMatrix.__name__,type(v).__name__))
