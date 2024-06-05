"""
Module defining basal square breeding value matrices and error checking routines.
"""

__all__ = [
    "SquareBreedingValueMatrix",
    "check_is_SquareBreedingValueMatrix",
]

from abc import ABCMeta, abstractmethod

import numpy
from pybrops.core.io.CSVInputOutput import CSVInputOutput
from pybrops.core.io.PandasInputOutput import PandasInputOutput
from pybrops.core.mat.ScaledSquareTaxaTraitMatrix import ScaledSquareTaxaTraitMatrix

class SquareBreedingValueMatrix(
        ScaledSquareTaxaTraitMatrix,
        PandasInputOutput,
        CSVInputOutput,
        metaclass = ABCMeta,
    ):
    """
    The SquareBreedingValueMatrix class represents a Multivariate Breeding Value.

    Notes
    -----
    All elements within a BreedingValueMatrix are mean-centered and scaled to
    unit variance for each trait.

    .. math::
        BV = \\frac{X - \\mu}{\\sigma}

    Where:

    - :math:`BV` is the breeding value.
    - :math:`X` is the phenotype value.
    - :math:`\\mu` is the mean (location) for :math:`X`.
    - :math:`\\sigma` is the standard deviation (scale) for :math:`X`.

    Phenotype values can be reconstituted using:

    .. math::
        X = \\sigma BV + \\mu
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################

    ########################## Special Object Methods ##########################

    ############## Forward numeric operators ###############
    ### __add__                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __sub__                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __mul__                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __matmul__              inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __truediv__             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __floordiv__            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __mod__                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __divmod__              inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __pow__                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __lshift__              inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rshift__              inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __and__                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __xor__                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __or__                  inherited from ``ScaledSquareTaxaTraitMatrix``

    ############# Backwards numeric operators ##############
    ### __radd__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rsub__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rmul__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rmatmul__             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rtruediv__            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rfloordiv__           inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rmod__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rdivmod__             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rlshift__             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rrshift__             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rand__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __rxor__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __ror__                 inherited from ``ScaledSquareTaxaTraitMatrix``

    ############# Augmented numeric operators ##############
    ### __iadd__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __isub__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __imul__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __imatmul__             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __itruediv__            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __ifloordiv__           inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __imod__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __ipow__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __ilshift__             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __irshift__             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __iand__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __ixor__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __ior__                 inherited from ``ScaledSquareTaxaTraitMatrix``

    ################## Logical operators ###################
    ### __lt__                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __le__                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __eq__                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __ne__                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __gt__                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __ge__                  inherited from ``ScaledSquareTaxaTraitMatrix``

    ################# Container operators ##################
    ### __len__                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __getitem__             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __setitem__             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __delitem__             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __iter__                inherited from ``ScaledSquareTaxaTraitMatrix``

    #################### Matrix copying ####################
    ### __copy__                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### __deepcopy__            inherited from ``ScaledSquareTaxaTraitMatrix``

    ########### Miscellaneous special functions ############
    ### __repr__                inherited from ``ScaledSquareTaxaTraitMatrix``

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    ### mat                     inherited from ``ScaledSquareTaxaTraitMatrix``

    ############## Matrix Metadata Properties ##############
    ### mat_ndim                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### mat_shape               inherited from ``ScaledSquareTaxaTraitMatrix``

    ############## Square Metadata Properties ##############
    ### nsquare                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### square_axes             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### square_axes_len         inherited from ``ScaledSquareTaxaTraitMatrix``

    ########### Square Taxa Metadata Properties ############
    ### nsquare_taxa            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### square_taxa_axes        inherited from ``ScaledSquareTaxaTraitMatrix``
    ### square_taxa_axes_len    inherited from ``ScaledSquareTaxaTraitMatrix``

    ################# Taxa Data Properites #################
    ### taxa                    inherited from ``ScaledSquareTaxaTraitMatrix``
    ### taxa_grp                inherited from ``ScaledSquareTaxaTraitMatrix``

    ############### Taxa Metadata Properites ###############
    ### ntaxa                   inherited from ``ScaledSquareTaxaTraitMatrix``
    ### taxa_axis               inherited from ``ScaledSquareTaxaTraitMatrix``
    ### taxa_grp_name           inherited from ``ScaledSquareTaxaTraitMatrix``
    ### taxa_grp_stix           inherited from ``ScaledSquareTaxaTraitMatrix``
    ### taxa_grp_spix           inherited from ``ScaledSquareTaxaTraitMatrix``
    ### taxa_grp_len            inherited from ``ScaledSquareTaxaTraitMatrix``

    ###################### Trait data ######################
    ### trait                   inherited from ``ScaledSquareTaxaTraitMatrix``

    #################### Trait metadata ####################
    ### ntrait                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### trait_axis              inherited from ``ScaledSquareTaxaTraitMatrix``

    ############## Scale metadata properties ###############
    ### location                inherited from ``ScaledSquareTaxaTraitMatrix``
    ### scale                   inherited from ``ScaledSquareTaxaTraitMatrix``

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    ### copy                    inherited from ``ScaledSquareTaxaTraitMatrix``
    ### deepcopy                inherited from ``ScaledSquareTaxaTraitMatrix``

    ######### Matrix element copy-on-manipulation ##########
    ### adjoin                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### delete                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### insert                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### select                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### adjoin_taxa             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### delete_taxa             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### insert_taxa             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### select_taxa             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### adjoin_trait            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### delete_trait            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### insert_trait            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### select_trait            inherited from ``ScaledSquareTaxaTraitMatrix``

    ######### Matrix element in-place-manipulation #########
    ### append                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### remove                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### incorp                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### append_taxa             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### remove_taxa             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### incorp_taxa             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### append_trait            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### remove_trait            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### incorp_trait            inherited from ``ScaledSquareTaxaTraitMatrix``

    ################### Sorting Methods ####################
    ### lexsort                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### reorder                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### sort                    inherited from ``ScaledSquareTaxaTraitMatrix``
    ### lexsort_taxa            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### reorder_taxa            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### sort_taxa               inherited from ``ScaledSquareTaxaTraitMatrix``
    ### lexsort_trait           inherited from ``ScaledSquareTaxaTraitMatrix``
    ### reorder_trait           inherited from ``ScaledSquareTaxaTraitMatrix``
    ### sort_trait              inherited from ``ScaledSquareTaxaTraitMatrix``

    ################### Grouping Methods ###################
    ### group                   inherited from ``ScaledSquareTaxaTraitMatrix``
    ### ungroup                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### is_grouped              inherited from ``ScaledSquareTaxaTraitMatrix``
    ### group_taxa              inherited from ``ScaledSquareTaxaTraitMatrix``
    ### ungroup_taxa            inherited from ``ScaledSquareTaxaTraitMatrix``
    ### is_grouped_taxa         inherited from ``ScaledSquareTaxaTraitMatrix``

    #################### Square Methods ####################
    ### is_square               inherited from ``ScaledSquareTaxaTraitMatrix``
    ### is_square_taxa          inherited from ``ScaledSquareTaxaTraitMatrix``

    ################### Scaling methods ####################
    ### transform               inherited from ``ScaledSquareTaxaTraitMatrix``
    ### untransform             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### rescale                 inherited from ``ScaledSquareTaxaTraitMatrix``
    ### unscale                 inherited from ``ScaledSquareTaxaTraitMatrix``

    ############################## Object Methods ##############################

    ############## Matrix summary statistics ###############
    @abstractmethod
    def targmax(self) -> numpy.ndarray:
        """
        Return indices of the maximum values along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing indices of maximum
            values along the trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def targmin(self) -> numpy.ndarray:
        """
        Return indices of the minimum values along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing indices of minimum
            values along the trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def tmax(self, unscale: bool) -> numpy.ndarray:
        """
        Return the maximum along the trait axis.

        Parameters
        ----------
        unscale : bool
            whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing maximum values along the trait
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def tmean(self, unscale: bool) -> numpy.ndarray:
        """
        Return the mean along the trait axis.

        Parameters
        ----------
        unscale : bool
            whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing maximum values along the trait
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def tmin(self, unscale: bool) -> numpy.ndarray:
        """
        Return the minimum along the trait axis.

        Parameters
        ----------
        unscale : bool
            whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing minimum values along the
            trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def trange(self, unscale: bool) -> numpy.ndarray:
        """
        Return the range along the trait axis.

        Parameters
        ----------
        unscale : bool
            whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing variance values along the
            trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def tstd(self, unscale: bool) -> numpy.ndarray:
        """
        Return the standard deviation along the trait axis.

        Parameters
        ----------
        unscale : bool
            whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing standard deviation values
            along the trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def tvar(self, unscale: bool) -> numpy.ndarray:
        """
        Return the variance along the trait axis.

        Parameters
        ----------
        unscale : bool
            whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing variance values along the
            trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    #################### Export Methods ####################
    ### to_pandas (inherited)   inherited from ``PandasInputOutput``
    ### to_csv (inherited)      inherited from ``CSVInputOutput``
    ### to_hdf5 (inherited)     inherited from ``ScaledSquareTaxaTraitMatrix``

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat                  inherited from ``ScaledSquareTaxaTraitMatrix``
    ### concat_taxa             inherited from ``ScaledSquareTaxaTraitMatrix``
    ### concat_trait            inherited from ``ScaledSquareTaxaTraitMatrix``

    #################### Import Methods ####################
    ### from_pandas             inherited from ``PandasInputOutput``
    ### from_csv                inherited from ``CSVInputOutput``
    ### from_hdf5               inherited from ``ScaledSquareTaxaTraitMatrix``



def check_is_SquareBreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``SquareBreedingValueMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, SquareBreedingValueMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                SquareBreedingValueMatrix.__name__,
                type(v).__name__
            )
        )
