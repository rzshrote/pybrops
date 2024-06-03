"""
Module defining interfaces and associated error checking routines for matrices
with taxa axes that are square and a trait axis which is not square.
"""

__all__ = [
    "SquareTaxaTraitMatrix",
    "check_is_SquareTaxaTraitMatrix",
]

from abc import ABCMeta
from pybrops.core.mat.SquareTaxaMatrix import SquareTaxaMatrix
from pybrops.core.mat.TraitMatrix import TraitMatrix


class SquareTaxaTraitMatrix(
        SquareTaxaMatrix,
        TraitMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for matrix wrapper object with taxa axes which are square
    and a trait axis which is not square.

    The purpose of this abstract class is to merge the following interfaces:

        1. SquareTaxaMatrix
        2. TraitMatrix
    """
    ########################## Special Object Methods ##########################

    ############## Forward numeric operators ###############
    ### __add__                 inherited from ``SquareTaxaMatrix``
    ### __sub__                 inherited from ``SquareTaxaMatrix``
    ### __mul__                 inherited from ``SquareTaxaMatrix``
    ### __matmul__              inherited from ``SquareTaxaMatrix``
    ### __truediv__             inherited from ``SquareTaxaMatrix``
    ### __floordiv__            inherited from ``SquareTaxaMatrix``
    ### __mod__                 inherited from ``SquareTaxaMatrix``
    ### __divmod__              inherited from ``SquareTaxaMatrix``
    ### __pow__                 inherited from ``SquareTaxaMatrix``
    ### __lshift__              inherited from ``SquareTaxaMatrix``
    ### __rshift__              inherited from ``SquareTaxaMatrix``
    ### __and__                 inherited from ``SquareTaxaMatrix``
    ### __xor__                 inherited from ``SquareTaxaMatrix``
    ### __or__                  inherited from ``SquareTaxaMatrix``

    ############# Backwards numeric operators ##############
    ### __radd__                inherited from ``SquareTaxaMatrix``
    ### __rsub__                inherited from ``SquareTaxaMatrix``
    ### __rmul__                inherited from ``SquareTaxaMatrix``
    ### __rmatmul__             inherited from ``SquareTaxaMatrix``
    ### __rtruediv__            inherited from ``SquareTaxaMatrix``
    ### __rfloordiv__           inherited from ``SquareTaxaMatrix``
    ### __rmod__                inherited from ``SquareTaxaMatrix``
    ### __rdivmod__             inherited from ``SquareTaxaMatrix``
    ### __rlshift__             inherited from ``SquareTaxaMatrix``
    ### __rrshift__             inherited from ``SquareTaxaMatrix``
    ### __rand__                inherited from ``SquareTaxaMatrix``
    ### __rxor__                inherited from ``SquareTaxaMatrix``
    ### __ror__                 inherited from ``SquareTaxaMatrix``

    ############# Augmented numeric operators ##############
    ### __iadd__                inherited from ``SquareTaxaMatrix``
    ### __isub__                inherited from ``SquareTaxaMatrix``
    ### __imul__                inherited from ``SquareTaxaMatrix``
    ### __imatmul__             inherited from ``SquareTaxaMatrix``
    ### __itruediv__            inherited from ``SquareTaxaMatrix``
    ### __ifloordiv__           inherited from ``SquareTaxaMatrix``
    ### __imod__                inherited from ``SquareTaxaMatrix``
    ### __ipow__                inherited from ``SquareTaxaMatrix``
    ### __ilshift__             inherited from ``SquareTaxaMatrix``
    ### __irshift__             inherited from ``SquareTaxaMatrix``
    ### __iand__                inherited from ``SquareTaxaMatrix``
    ### __ixor__                inherited from ``SquareTaxaMatrix``
    ### __ior__                 inherited from ``SquareTaxaMatrix``

    ################## Logical operators ###################
    ### __lt__                  inherited from ``SquareTaxaMatrix``
    ### __le__                  inherited from ``SquareTaxaMatrix``
    ### __eq__                  inherited from ``SquareTaxaMatrix``
    ### __ne__                  inherited from ``SquareTaxaMatrix``
    ### __gt__                  inherited from ``SquareTaxaMatrix``
    ### __ge__                  inherited from ``SquareTaxaMatrix``

    ################# Container operators ##################
    ### __len__                 inherited from ``SquareTaxaMatrix``
    ### __getitem__             inherited from ``SquareTaxaMatrix``
    ### __setitem__             inherited from ``SquareTaxaMatrix``
    ### __delitem__             inherited from ``SquareTaxaMatrix``
    ### __iter__                inherited from ``SquareTaxaMatrix``

    #################### Matrix copying ####################
    ### __copy__                inherited from ``SquareTaxaMatrix``
    ### __deepcopy__            inherited from ``SquareTaxaMatrix``

    ########### Miscellaneous special functions ############
    ### __repr__                inherited from ``SquareTaxaMatrix``

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    ### mat                     inherited from ``SquareTaxaMatrix``

    ############## Matrix Metadata Properties ##############
    ### mat_ndim                inherited from ``SquareTaxaMatrix``
    ### mat_shape               inherited from ``SquareTaxaMatrix``

    ############## Square Metadata Properties ##############
    ### nsquare                 inherited from ``SquareTaxaMatrix``
    ### square_axes             inherited from ``SquareTaxaMatrix``
    ### square_axes_len         inherited from ``SquareTaxaMatrix``

    ########### Square Taxa Metadata Properties ############
    ### nsquare_taxa            inherited from ``SquareTaxaMatrix``
    ### square_taxa_axes        inherited from ``SquareTaxaMatrix``
    ### square_taxa_axes_len    inherited from ``SquareTaxaMatrix``

    ################# Taxa Data Properites #################
    ### taxa                    inherited from ``SquareTaxaMatrix``
    ### taxa_grp                inherited from ``SquareTaxaMatrix``

    ############### Taxa Metadata Properites ###############
    ### ntaxa                   inherited from ``SquareTaxaMatrix``
    ### taxa_axis               inherited from ``SquareTaxaMatrix``
    ### taxa_grp_name           inherited from ``SquareTaxaMatrix``
    ### taxa_grp_stix           inherited from ``SquareTaxaMatrix``
    ### taxa_grp_spix           inherited from ``SquareTaxaMatrix``
    ### taxa_grp_len            inherited from ``SquareTaxaMatrix``

    ###################### Trait data ######################
    ### trait                   inherited from ``TraitMatrix``

    #################### Trait metadata ####################
    ### ntrait                  inherited from ``TraitMatrix``
    ### trait_axis              inherited from ``TraitMatrix``

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    ### copy                    inherited from ``SquareTaxaMatrix``
    ### deepcopy                inherited from ``SquareTaxaMatrix``

    ######### Matrix element copy-on-manipulation ##########
    ### adjoin                  inherited from ``SquareTaxaMatrix``
    ### delete                  inherited from ``SquareTaxaMatrix``
    ### insert                  inherited from ``SquareTaxaMatrix``
    ### select                  inherited from ``SquareTaxaMatrix``
    ### adjoin_taxa             inherited from ``SquareTaxaMatrix``
    ### delete_taxa             inherited from ``SquareTaxaMatrix``
    ### insert_taxa             inherited from ``SquareTaxaMatrix``
    ### select_taxa             inherited from ``SquareTaxaMatrix``
    ### adjoin_trait            inherited from ``TraitMatrix``
    ### delete_trait            inherited from ``TraitMatrix``
    ### insert_trait            inherited from ``TraitMatrix``
    ### select_trait            inherited from ``TraitMatrix``

    ######### Matrix element in-place-manipulation #########
    ### append                  inherited from ``SquareTaxaMatrix``
    ### remove                  inherited from ``SquareTaxaMatrix``
    ### incorp                  inherited from ``SquareTaxaMatrix``
    ### append_taxa             inherited from ``SquareTaxaMatrix``
    ### remove_taxa             inherited from ``SquareTaxaMatrix``
    ### incorp_taxa             inherited from ``SquareTaxaMatrix``
    ### append_trait            inherited from ``TraitMatrix``
    ### remove_trait            inherited from ``TraitMatrix``
    ### incorp_trait            inherited from ``TraitMatrix``

    ################### Sorting Methods ####################
    ### lexsort                 inherited from ``SquareTaxaMatrix``
    ### reorder                 inherited from ``SquareTaxaMatrix``
    ### sort                    inherited from ``SquareTaxaMatrix``
    ### lexsort_taxa            inherited from ``SquareTaxaMatrix``
    ### reorder_taxa            inherited from ``SquareTaxaMatrix``
    ### sort_taxa               inherited from ``SquareTaxaMatrix``
    ### lexsort_trait           inherited from ``TraitMatrix``
    ### reorder_trait           inherited from ``TraitMatrix``
    ### sort_trait              inherited from ``TraitMatrix``

    ################### Grouping Methods ###################
    ### group                   inherited from ``SquareTaxaMatrix``
    ### ungroup                 inherited from ``SquareTaxaMatrix``
    ### is_grouped              inherited from ``SquareTaxaMatrix``
    ### group_taxa              inherited from ``SquareTaxaMatrix``
    ### ungroup_taxa            inherited from ``SquareTaxaMatrix``
    ### is_grouped_taxa         inherited from ``SquareTaxaMatrix``

    #################### Square Methods ####################
    ### is_square               inherited from ``SquareTaxaMatrix``
    ### is_square_taxa          inherited from ``SquareTaxaMatrix``

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat                  inherited from ``SquareTaxaMatrix``
    ### concat_taxa             inherited from ``SquareTaxaMatrix``
    ### concat_trait            inherited from ``TraitMatrix``
    pass



################################## Utilities ###################################
def check_is_SquareTaxaTraitMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``SquareTaxaTraitMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, SquareTaxaTraitMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,SquareTaxaTraitMatrix.__name__,type(v).__name__))
