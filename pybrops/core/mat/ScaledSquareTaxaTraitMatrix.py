"""
Module defining interfaces and associated error checking routines for scaled 
matrices with taxa axes that are square and a trait axis which is not square.
"""


from abc import ABCMeta
from pybrops.core.mat.ScaledMatrix import ScaledMatrix
from pybrops.core.mat.SquareTaxaTraitMatrix import SquareTaxaTraitMatrix


class ScaledSquareTaxaTraitMatrix(
        SquareTaxaTraitMatrix,
        ScaledMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for scaled matrix wrapper object with taxa axes which are
    square and a trait axis which is not square.

    The purpose of this abstract class is to merge the following interfaces:

        1. SquareTaxaTraitMatrix
        2. ScaledMatrix
    """
    ########################## Special Object Methods ##########################

    ############## Forward numeric operators ###############
    ### __add__                 inherited from ``SquareTaxaTraitMatrix``
    ### __sub__                 inherited from ``SquareTaxaTraitMatrix``
    ### __mul__                 inherited from ``SquareTaxaTraitMatrix``
    ### __matmul__              inherited from ``SquareTaxaTraitMatrix``
    ### __truediv__             inherited from ``SquareTaxaTraitMatrix``
    ### __floordiv__            inherited from ``SquareTaxaTraitMatrix``
    ### __mod__                 inherited from ``SquareTaxaTraitMatrix``
    ### __divmod__              inherited from ``SquareTaxaTraitMatrix``
    ### __pow__                 inherited from ``SquareTaxaTraitMatrix``
    ### __lshift__              inherited from ``SquareTaxaTraitMatrix``
    ### __rshift__              inherited from ``SquareTaxaTraitMatrix``
    ### __and__                 inherited from ``SquareTaxaTraitMatrix``
    ### __xor__                 inherited from ``SquareTaxaTraitMatrix``
    ### __or__                  inherited from ``SquareTaxaTraitMatrix``

    ############# Backwards numeric operators ##############
    ### __radd__                inherited from ``SquareTaxaTraitMatrix``
    ### __rsub__                inherited from ``SquareTaxaTraitMatrix``
    ### __rmul__                inherited from ``SquareTaxaTraitMatrix``
    ### __rmatmul__             inherited from ``SquareTaxaTraitMatrix``
    ### __rtruediv__            inherited from ``SquareTaxaTraitMatrix``
    ### __rfloordiv__           inherited from ``SquareTaxaTraitMatrix``
    ### __rmod__                inherited from ``SquareTaxaTraitMatrix``
    ### __rdivmod__             inherited from ``SquareTaxaTraitMatrix``
    ### __rlshift__             inherited from ``SquareTaxaTraitMatrix``
    ### __rrshift__             inherited from ``SquareTaxaTraitMatrix``
    ### __rand__                inherited from ``SquareTaxaTraitMatrix``
    ### __rxor__                inherited from ``SquareTaxaTraitMatrix``
    ### __ror__                 inherited from ``SquareTaxaTraitMatrix``

    ############# Augmented numeric operators ##############
    ### __iadd__                inherited from ``SquareTaxaTraitMatrix``
    ### __isub__                inherited from ``SquareTaxaTraitMatrix``
    ### __imul__                inherited from ``SquareTaxaTraitMatrix``
    ### __imatmul__             inherited from ``SquareTaxaTraitMatrix``
    ### __itruediv__            inherited from ``SquareTaxaTraitMatrix``
    ### __ifloordiv__           inherited from ``SquareTaxaTraitMatrix``
    ### __imod__                inherited from ``SquareTaxaTraitMatrix``
    ### __ipow__                inherited from ``SquareTaxaTraitMatrix``
    ### __ilshift__             inherited from ``SquareTaxaTraitMatrix``
    ### __irshift__             inherited from ``SquareTaxaTraitMatrix``
    ### __iand__                inherited from ``SquareTaxaTraitMatrix``
    ### __ixor__                inherited from ``SquareTaxaTraitMatrix``
    ### __ior__                 inherited from ``SquareTaxaTraitMatrix``

    ################## Logical operators ###################
    ### __lt__                  inherited from ``SquareTaxaTraitMatrix``
    ### __le__                  inherited from ``SquareTaxaTraitMatrix``
    ### __eq__                  inherited from ``SquareTaxaTraitMatrix``
    ### __ne__                  inherited from ``SquareTaxaTraitMatrix``
    ### __gt__                  inherited from ``SquareTaxaTraitMatrix``
    ### __ge__                  inherited from ``SquareTaxaTraitMatrix``

    ################# Container operators ##################
    ### __len__                 inherited from ``SquareTaxaTraitMatrix``
    ### __getitem__             inherited from ``SquareTaxaTraitMatrix``
    ### __setitem__             inherited from ``SquareTaxaTraitMatrix``
    ### __delitem__             inherited from ``SquareTaxaTraitMatrix``
    ### __iter__                inherited from ``SquareTaxaTraitMatrix``

    #################### Matrix copying ####################
    ### __copy__                inherited from ``SquareTaxaTraitMatrix``
    ### __deepcopy__            inherited from ``SquareTaxaTraitMatrix``

    ########### Miscellaneous special functions ############
    ### __repr__                inherited from ``SquareTaxaTraitMatrix``

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    ### mat                     inherited from ``SquareTaxaTraitMatrix``

    ############## Matrix Metadata Properties ##############
    ### mat_ndim                inherited from ``SquareTaxaTraitMatrix``
    ### mat_shape               inherited from ``SquareTaxaTraitMatrix``

    ############## Square Metadata Properties ##############
    ### nsquare                 inherited from ``SquareTaxaTraitMatrix``
    ### square_axes             inherited from ``SquareTaxaTraitMatrix``
    ### square_axes_len         inherited from ``SquareTaxaTraitMatrix``

    ########### Square Taxa Metadata Properties ############
    ### nsquare_taxa            inherited from ``SquareTaxaTraitMatrix``
    ### square_taxa_axes        inherited from ``SquareTaxaTraitMatrix``
    ### square_taxa_axes_len    inherited from ``SquareTaxaTraitMatrix``

    ################# Taxa Data Properites #################
    ### taxa                    inherited from ``SquareTaxaTraitMatrix``
    ### taxa_grp                inherited from ``SquareTaxaTraitMatrix``

    ############### Taxa Metadata Properites ###############
    ### ntaxa                   inherited from ``SquareTaxaTraitMatrix``
    ### taxa_axis               inherited from ``SquareTaxaTraitMatrix``
    ### taxa_grp_name           inherited from ``SquareTaxaTraitMatrix``
    ### taxa_grp_stix           inherited from ``SquareTaxaTraitMatrix``
    ### taxa_grp_spix           inherited from ``SquareTaxaTraitMatrix``
    ### taxa_grp_len            inherited from ``SquareTaxaTraitMatrix``

    ###################### Trait data ######################
    ### trait                   inherited from ``SquareTaxaTraitMatrix``

    #################### Trait metadata ####################
    ### ntrait                  inherited from ``SquareTaxaTraitMatrix``
    ### trait_axis              inherited from ``SquareTaxaTraitMatrix``

    ############## Scale metadata properties ###############
    ### location                inherited from ``ScaledMatrix``
    ### scale                   inherited from ``ScaledMatrix``

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    ### copy                    inherited from ``SquareTaxaTraitMatrix``
    ### deepcopy                inherited from ``SquareTaxaTraitMatrix``

    ######### Matrix element copy-on-manipulation ##########
    ### adjoin                  inherited from ``SquareTaxaTraitMatrix``
    ### delete                  inherited from ``SquareTaxaTraitMatrix``
    ### insert                  inherited from ``SquareTaxaTraitMatrix``
    ### select                  inherited from ``SquareTaxaTraitMatrix``
    ### adjoin_taxa             inherited from ``SquareTaxaTraitMatrix``
    ### delete_taxa             inherited from ``SquareTaxaTraitMatrix``
    ### insert_taxa             inherited from ``SquareTaxaTraitMatrix``
    ### select_taxa             inherited from ``SquareTaxaTraitMatrix``
    ### adjoin_trait            inherited from ``SquareTaxaTraitMatrix``
    ### delete_trait            inherited from ``SquareTaxaTraitMatrix``
    ### insert_trait            inherited from ``SquareTaxaTraitMatrix``
    ### select_trait            inherited from ``SquareTaxaTraitMatrix``

    ######### Matrix element in-place-manipulation #########
    ### append                  inherited from ``SquareTaxaTraitMatrix``
    ### remove                  inherited from ``SquareTaxaTraitMatrix``
    ### incorp                  inherited from ``SquareTaxaTraitMatrix``
    ### append_taxa             inherited from ``SquareTaxaTraitMatrix``
    ### remove_taxa             inherited from ``SquareTaxaTraitMatrix``
    ### incorp_taxa             inherited from ``SquareTaxaTraitMatrix``
    ### append_trait            inherited from ``SquareTaxaTraitMatrix``
    ### remove_trait            inherited from ``SquareTaxaTraitMatrix``
    ### incorp_trait            inherited from ``SquareTaxaTraitMatrix``

    ################### Sorting Methods ####################
    ### lexsort                 inherited from ``SquareTaxaTraitMatrix``
    ### reorder                 inherited from ``SquareTaxaTraitMatrix``
    ### sort                    inherited from ``SquareTaxaTraitMatrix``
    ### lexsort_taxa            inherited from ``SquareTaxaTraitMatrix``
    ### reorder_taxa            inherited from ``SquareTaxaTraitMatrix``
    ### sort_taxa               inherited from ``SquareTaxaTraitMatrix``
    ### lexsort_trait           inherited from ``SquareTaxaTraitMatrix``
    ### reorder_trait           inherited from ``SquareTaxaTraitMatrix``
    ### sort_trait              inherited from ``SquareTaxaTraitMatrix``

    ################### Grouping Methods ###################
    ### group                   inherited from ``SquareTaxaTraitMatrix``
    ### ungroup                 inherited from ``SquareTaxaTraitMatrix``
    ### is_grouped              inherited from ``SquareTaxaTraitMatrix``
    ### group_taxa              inherited from ``SquareTaxaTraitMatrix``
    ### ungroup_taxa            inherited from ``SquareTaxaTraitMatrix``
    ### is_grouped_taxa         inherited from ``SquareTaxaTraitMatrix``

    #################### Square Methods ####################
    ### is_square               inherited from ``SquareTaxaTraitMatrix``
    ### is_square_taxa          inherited from ``SquareTaxaTraitMatrix``

    ################### Scaling methods ####################
    ### transform               inherited from ``ScaledMatrix``
    ### untransform             inherited from ``ScaledMatrix``
    ### rescale                 inherited from ``ScaledMatrix``
    ### unscale                 inherited from ``ScaledMatrix``

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat                  inherited from ``SquareTaxaTraitMatrix``
    ### concat_taxa             inherited from ``SquareTaxaTraitMatrix``
    ### concat_trait            inherited from ``SquareTaxaTraitMatrix``
    pass



################################## Utilities ###################################
def check_is_ScaledSquareTaxaTraitMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``ScaledSquareTaxaTraitMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, ScaledSquareTaxaTraitMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                ScaledSquareTaxaTraitMatrix.__name__,
                type(v).__name__
            )
        )
