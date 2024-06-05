"""
Module implementing dense square breeding value matrices.
"""

__all__ = [
    "DenseSquareBreedingValueMatrix",
    "check_is_DenseSquareBreedingValueMatrix",
]

import copy
from functools import reduce
from numbers import Integral, Real
from typing import Optional, Sequence, Union
import numpy
import pandas
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_type_python import check_Sequence_all_type, check_is_Integral, check_is_bool, check_is_str
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_column_indices
from pybrops.core.error.error_value_python import check_len
from pybrops.core.mat.DenseScaledSquareTaxaTraitMatrix import DenseScaledSquareTaxaTraitMatrix
from pybrops.core.util.array import flattenix
from pybrops.popgen.sbvmat import SquareBreedingValueMatrix


class DenseSquareBreedingValueMatrix(
        DenseScaledSquareTaxaTraitMatrix,
        SquareBreedingValueMatrix,
    ):
    """
    docstring for DenseSquareBreedingValueMatrix.
    """

    ########################## Special Object Methods ##########################

    ##################### Constructor ######################
    def __init__(
            self,
            mat: numpy.ndarray, 
            location: Union[numpy.ndarray,Real] = 0.0, 
            scale: Union[numpy.ndarray,Real] = 1.0, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseSquareBreedingValueMatrix.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # call DenseScaledSquareTaxaTraitMatrix constructor
        super(DenseSquareBreedingValueMatrix, self).__init__(
            mat = mat,
            location = location,
            scale = scale,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    ############## Forward numeric operators ###############
    ### __add__                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __sub__                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __mul__                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __matmul__              inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __truediv__             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __floordiv__            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __mod__                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __divmod__              inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __pow__                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __lshift__              inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rshift__              inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __and__                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __xor__                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __or__                  inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ############# Backwards numeric operators ##############
    ### __radd__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rsub__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rmul__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rmatmul__             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rtruediv__            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rfloordiv__           inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rmod__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rdivmod__             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rlshift__             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rrshift__             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rand__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __rxor__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __ror__                 inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ############# Augmented numeric operators ##############
    ### __iadd__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __isub__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __imul__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __imatmul__             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __itruediv__            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __ifloordiv__           inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __imod__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __ipow__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __ilshift__             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __irshift__             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __iand__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __ixor__                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __ior__                 inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ################## Logical operators ###################
    ### __lt__                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __le__                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __eq__                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __ne__                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __gt__                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __ge__                  inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ################# Container operators ##################
    ### __len__                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __getitem__             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __setitem__             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __delitem__             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### __iter__                inherited from ``DenseScaledSquareTaxaTraitMatrix``

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseSquareBreedingValueMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseSquareBreedingValueMatrix
            A copy of the DenseSquareBreedingValueMatrix.
        """
        return self.__class__(
            mat = copy.copy(self.mat),
            location = copy.copy(self.location),
            scale = copy.copy(self.scale),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
            trait = copy.copy(self.trait),
        )
    
    def __deepcopy__(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseSquareBreedingValueMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict, None, default = None
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseSquareBreedingValueMatrix
            A deep copy of the DenseSquareBreedingValueMatrix.
        """
        return self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            location = copy.deepcopy(self.location, memo),
            scale = copy.deepcopy(self.scale, memo),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo),
            trait = copy.deepcopy(self.trait, memo),
        )

    ########### Miscellaneous special functions ############
    ### __repr__                inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    ### mat                     inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ############## Matrix Metadata Properties ##############
    ### mat_ndim                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### mat_shape               inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ############## Square Metadata Properties ##############
    ### nsquare                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### square_axes             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### square_axes_len         inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ########### Square Taxa Metadata Properties ############
    ### nsquare_taxa            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### square_taxa_axes        inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### square_taxa_axes_len    inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ################# Taxa Data Properites #################
    ### taxa                    inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### taxa_grp                inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ############### Taxa Metadata Properites ###############
    ### ntaxa                   inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### taxa_axis               inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### taxa_grp_name           inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### taxa_grp_stix           inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### taxa_grp_spix           inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### taxa_grp_len            inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ###################### Trait data ######################
    ### trait                   inherited from ``DenseScaledSquareTaxaTraitMatrix``

    #################### Trait metadata ####################
    ### ntrait                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### trait_axis              inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ############## Scale metadata properties ###############
    ### location                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### scale                   inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseSquareBreedingValueMatrix':
        """
        Make a shallow copy of the DenseSquareBreedingValueMatrix.

        Returns
        -------
        out : DenseSquareBreedingValueMatrix
            A shallow copy of the original DenseSquareBreedingValueMatrix.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseSquareBreedingValueMatrix':
        """
        Make a deep copy of the DenseSquareBreedingValueMatrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseSquareBreedingValueMatrix
            A deep copy of the original DenseSquareBreedingValueMatrix.
        """
        return copy.deepcopy(self, memo)

    ######### Matrix element copy-on-manipulation ##########
    ### adjoin                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### delete                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### insert                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### select                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### adjoin_taxa             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### delete_taxa             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### insert_taxa             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### select_taxa             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### adjoin_trait            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### delete_trait            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### insert_trait            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### select_trait            inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ######### Matrix element in-place-manipulation #########
    ### append                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### remove                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### incorp                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### append_taxa             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### remove_taxa             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### incorp_taxa             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### append_trait            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### remove_trait            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### incorp_trait            inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ################### Sorting Methods ####################
    ### lexsort                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### reorder                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### sort                    inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### lexsort_taxa            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### reorder_taxa            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### sort_taxa               inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### lexsort_trait           inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### reorder_trait           inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### sort_trait              inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ################### Grouping Methods ###################
    ### group                   inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### ungroup                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### is_grouped              inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### group_taxa              inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### ungroup_taxa            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### is_grouped_taxa         inherited from ``DenseScaledSquareTaxaTraitMatrix``

    #################### Square Methods ####################
    ### is_square               inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### is_square_taxa          inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ################### Scaling methods ####################
    ### transform               inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### untransform             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### rescale                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### unscale                 inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ############################## Object Methods ##############################

    ############## Matrix summary statistics ###############
    def targmax(self) -> numpy.ndarray:
        """
        Return raveled indices of the maximum values along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing indices of maximum
            values along the trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        # get number of traits, dimensions, trait axis
        ntrait = self.ntrait
        ndim = self.mat_ndim
        trait_axis = self.trait_axis

        # get matrix
        mat = self.mat

        # allocate memory for empty array
        # (t,)
        out = numpy.empty((ntrait,), dtype=int)

        # for each trait
        for trait in range(ntrait):
            # construct a view index tuple
            ix = tuple(slice(None) if i!=trait_axis else trait for i in range(ndim))
            # get raveled index of maximum value
            out[trait] = mat[ix].argmax()

        return out

    def targmin(self) -> numpy.ndarray:
        """
        Return raveled indices of the minimum values along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing indices of minimum
            values along the trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        # get number of traits, dimensions, trait axis
        ntrait = self.ntrait
        ndim = self.mat_ndim
        trait_axis = self.trait_axis

        # get matrix
        mat = self.mat

        # allocate memory for empty array
        # (t,)
        out = numpy.empty((ntrait,), dtype=int)

        # for each trait
        for trait in range(ntrait):
            # construct a view index tuple
            ix = tuple(slice(None) if i!=trait_axis else trait for i in range(ndim))
            # get raveled index of maximum value
            out[trait] = mat[ix].argmin()

        return out

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
        # type checks
        check_is_bool(unscale, "unscale")

        # get number of traits, dimensions, trait axis
        ntrait = self.ntrait
        ndim = self.mat_ndim
        trait_axis = self.trait_axis

        # get matrix
        mat = self.mat

        # allocate memory for empty array
        # (t,)
        out = numpy.empty((ntrait,), dtype=mat.dtype)

        # for each trait
        for trait in range(ntrait):
            # construct a view index tuple
            ix = tuple(slice(None) if i!=trait_axis else trait for i in range(ndim))
            # get raveled maximum value
            out[trait] = mat[ix].max()

        # unscale results if needed
        # (t,) *= (t,)
        # (t,) += (t,)
        if unscale:
            out *= self.scale
            out += self.location

        return out

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
        # type checks
        check_is_bool(unscale, "unscale")

        # get number of traits, dimensions, trait axis
        ntrait = self.ntrait
        ndim = self.mat_ndim
        trait_axis = self.trait_axis

        # get matrix
        mat = self.mat

        # allocate memory for empty array
        # (t,)
        out = numpy.empty((ntrait,), dtype=mat.dtype)

        # for each trait
        for trait in range(ntrait):
            # construct a view index tuple
            ix = tuple(slice(None) if i!=trait_axis else trait for i in range(ndim))
            # get raveled mean value
            out[trait] = mat[ix].mean()

        # unscale results if needed
        # (t,) *= (t,)
        # (t,) += (t,)
        if unscale:
            out *= self.scale
            out += self.location

        return out

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
        # type checks
        check_is_bool(unscale, "unscale")

        # get number of traits, dimensions, trait axis
        ntrait = self.ntrait
        ndim = self.mat_ndim
        trait_axis = self.trait_axis

        # get matrix
        mat = self.mat

        # allocate memory for empty array
        # (t,)
        out = numpy.empty((ntrait,), dtype=mat.dtype)

        # for each trait
        for trait in range(ntrait):
            # construct a view index tuple
            ix = tuple(slice(None) if i!=trait_axis else trait for i in range(ndim))
            # get raveled minimum value
            out[trait] = mat[ix].min()

        # unscale results if needed
        # (t,) *= (t,)
        # (t,) += (t,)
        if unscale:
            out *= self.scale
            out += self.location

        return out

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
        # type checks
        check_is_bool(unscale, "unscale")

        # get number of traits, dimensions, trait axis
        ntrait = self.ntrait
        ndim = self.mat_ndim
        trait_axis = self.trait_axis

        # get matrix
        mat = self.mat

        # allocate memory for empty array
        # (t,)
        out = numpy.empty((ntrait,), dtype=mat.dtype)

        # for each trait
        for trait in range(ntrait):
            # construct a view index tuple
            ix = tuple(slice(None) if i!=trait_axis else trait for i in range(ndim))
            # get raveled minimum value
            out[trait] = mat[ix].ptp()

        # unscale results if needed
        # (t,) *= (t,)
        if unscale:
            out *= self.scale

        return out

    # TODO: use unscale argument
    def tstd(self, unscale: bool) -> numpy.ndarray:
        """
        Return the standard deviation along the trait axis.

        Parameters
        ----------
        unscale : bool
            Unused argument.
            whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing standard deviation values
            along the trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        # type checks
        check_is_bool(unscale, "unscale")

        # get number of traits, dimensions, trait axis
        ntrait = self.ntrait
        ndim = self.mat_ndim
        trait_axis = self.trait_axis

        # get unscaled copy of matrix
        mat = self.unscale(inplace = False)
        # allocate memory for empty array
        # (t,)
        out = numpy.empty((ntrait,), dtype=mat.dtype)

        # for each trait
        for trait in range(ntrait):
            # construct a view index tuple
            ix = tuple(slice(None) if i!=trait_axis else trait for i in range(ndim))
            # get raveled standard deviation value
            out[trait] = mat[ix].std()

        return out

    # TODO: use unscale argument
    def tvar(self, unscale: bool) -> numpy.ndarray:
        """
        Return the variance along the trait axis.

        Parameters
        ----------
        unscale : bool
            Unused argument.
            whether to transform results to their unscaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing variance values along the
            trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        # type checks
        check_is_bool(unscale, "unscale")

        # get number of traits, dimensions, trait axis
        ntrait = self.ntrait
        ndim = self.mat_ndim
        trait_axis = self.trait_axis

        # get unscaled copy of matrix
        mat = self.unscale(inplace = False)
        # allocate memory for empty array
        # (t,)
        out = numpy.empty((ntrait,), dtype=mat.dtype)

        # for each trait
        for trait in range(ntrait):
            # construct a view index tuple
            ix = tuple(slice(None) if i!=trait_axis else trait for i in range(ndim))
            # get raveled standard deviation value
            out[trait] = mat[ix].var()

        return out

    #################### Export Methods ####################
    ### to_hdf5 (inherited)     inherited from ``DenseScaledSquareTaxaTraitMatrix``

    def to_pandas(
            self,
            taxa_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            taxa_grp_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            trait_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            value_colname: str = "value",
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Export a DenseSquareBreedingValueMatrix to a pandas.DataFrame.

        Parameters
        ----------
        taxa_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name taxa axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one taxa axis column is assumed.
            If ``None``, then do not export taxa columns.
        
        taxa_grp_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name taxa group axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one taxa group axis column is assumed.
            If ``None``, then do not export taxa group columns.

        trait_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name trait columns.
            There must be one column name for each trait.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one column name is assumed.
            If ``None``, then do not export trait column.
        
        value_colname : str, default = "value"
            Name of the value column.
        
        kwargs : dict
            Additional keyword arguments to be passed for dataframe construction.

        Returns
        -------
        out : pandas.DataFrame
            Output pandas dataframe.
        """
        ###
        ### Process inputs
        ###

        # get the number of taxa axes
        ntaxaaxes = self.nsquare_taxa

        # get the number of trait axes
        ntraitaxes = 1

        ### process inputs: taxa_colnames

        # process ``None``
        if taxa_colnames is None:
            taxa_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(taxa_colnames, bool):
            taxa_colnames = ["taxa_"+str(i).zfill(len(str(ntaxaaxes))) for i in range(ntaxaaxes)] if taxa_colnames else [None for i in range(ntaxaaxes)]
        elif isinstance(taxa_colnames, str):
            taxa_colnames = [taxa_colnames]
        elif isinstance(taxa_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``taxa_colnames`` must be of type ``bool``, ``str``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(taxa_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(taxa_colnames, "taxa_colnames", ntaxaaxes)

        # check sequence element types
        check_Sequence_all_type(taxa_colnames, "taxa_colnames", (str,type(None)))
        
        ### process inputs: taxa_grp_colnames

        # process ``None``
        if taxa_grp_colnames is None:
            taxa_grp_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(taxa_grp_colnames, bool):
            taxa_grp_colnames = ["taxa_grp_"+str(i).zfill(len(str(ntaxaaxes))) for i in range(ntaxaaxes)] if taxa_grp_colnames else [None for i in range(ntaxaaxes)]
        elif isinstance(taxa_grp_colnames, str):
            taxa_grp_colnames = [taxa_grp_colnames]
        elif isinstance(taxa_grp_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``taxa_grp_colnames`` must be of type ``bool``, ``str``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(taxa_grp_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(taxa_grp_colnames, "taxa_grp_colnames", ntaxaaxes)

        # check sequence element types
        check_Sequence_all_type(taxa_grp_colnames, "taxa_grp_colnames", (str,type(None)))

        ### process inputs: trait_colnames

        # process ``None``
        if trait_colnames is None:
            trait_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(trait_colnames, bool):
            trait_colnames = ["trait_"+str(i).zfill(len(str(ntraitaxes))) for i in range(ntraitaxes)] if trait_colnames else [None for i in range(ntraitaxes)]
        elif isinstance(trait_colnames, str):
            trait_colnames = [trait_colnames]
        elif isinstance(trait_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``trait_colnames`` must be of type ``bool``, ``str``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(trait_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(trait_colnames, "trait_colnames", ntraitaxes)

        # check sequence element types
        check_Sequence_all_type(trait_colnames, "trait_colnames", (str,type(None)))

        ### process inputs: value_colname
        check_is_str(value_colname, "value_colname")

        ###
        ### Process data
        ###

        # get taxa names
        taxa = numpy.array(["Taxon"+str(e).zfill(len(str(self.ntaxa))) for e in range(self.ntaxa)], dtype = object) if self.taxa is None else self.taxa

        # get taxa_grp names
        taxa_grp = numpy.array([pandas.NA for i in range(self.ntaxa)], dtype=object) if self.taxa_grp is None else self.taxa_grp

        # get trait names
        trait = numpy.array(["Trait"+str(i).zfill(len(str(self.ntrait))) for i in range(self.ntrait)], dtype=object) if self.trait is None else self.trait

        print(self.ntrait)

        # calculate flattened array and corresponding axis indices
        flatmat, axisix_data = flattenix(self.mat)

        ###
        ### Export data to dict then to pandas and return
        ###

        # create empty dict
        out_dict = {}

        # save taxa column data
        for ix,taxa_colname in zip(self.square_taxa_axes,taxa_colnames):
            if taxa_colname is None:
                continue
            axisix = axisix_data[ix]
            out_dict.update({taxa_colname: taxa[axisix]})

        # save taxa grp column data
        for ix,taxa_grp_colname in zip(self.square_taxa_axes,taxa_grp_colnames):
            if taxa_grp_colname is None:
                continue
            axisix = axisix_data[ix]
            out_dict.update({taxa_grp_colname: taxa_grp[axisix]})
        
        # save trait data
        trait_colname = trait_colnames[0]
        if trait_colname is not None:
            axisix = axisix_data[self.trait_axis]
            out_dict.update({trait_colname: trait[axisix]})
        
        # save values
        out_dict.update({value_colname: flatmat})

        # create a pandas DataFrame from the data
        out = pandas.DataFrame(out_dict, **kwargs)

        return out

    def to_csv(
            self,
            filename: str,
            taxa_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            taxa_grp_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            trait_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            value_colname: str = "value",
            sep: str = ',', 
            header: bool = True, 
            index: bool = False, 
            **kwargs: dict
        ) -> None:
        """
        Export a DenseSquareBreedingValueMatrix to a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name to which to write.
        
        taxa_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name taxa axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one taxa axis column is assumed.
            If ``None``, then do not export taxa columns.
        
        taxa_grp_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name taxa group axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one taxa group axis column is assumed.
            If ``None``, then do not export taxa group columns.

        trait_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name trait columns.
            There must be one column name for each trait.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one column name is assumed.
            If ``None``, then do not export trait column.
        
        value_colname : str, default = "value"
            Name of the value column.
        
        sep : str, default = ","
            Separator to use in the exported CSV file.
        
        header : bool, default = True
            Whether to save header names.
        
        index : bool, default = False
            Whether to save a row index in the exported CSV file.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a CSV.
        """
        # convert DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix to pandas.DataFrame
        df = self.to_pandas(
            taxa_colnames = taxa_colnames,
            taxa_grp_colnames = taxa_grp_colnames,
            trait_colnames = trait_colnames,
            value_colname = value_colname,
        )

        # export using pandas
        df.to_csv(
            path_or_buf = filename,
            sep         = sep,
            header      = header,
            index       = index,
            **kwargs
        )

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### concat_taxa             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### concat_trait            inherited from ``DenseScaledSquareTaxaTraitMatrix``

    #################### Import Methods ####################
    ### from_hdf5               inherited from ``DenseScaledSquareTaxaTraitMatrix``

    @classmethod
    def from_pandas(
            cls,
            df: pandas.DataFrame,
            taxa_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            taxa_grp_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            trait_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            value_colname: Union[str,Integral] = "value",
            ntaxaaxes: Integral = 2,
            **kwargs: dict
        ) -> 'DenseSquareBreedingValueMatrix':
        """
        Import a DenseSquareBreedingValueMatrix from a pandas.DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.

        taxa_colnames : bool, str, Integral, Sequence of str or Integral, None, default = True
            Sequence of column names or indices from which to read taxa axis columns.
            There must be one column name or index for each taxa axis.
            If an element in the sequence cannot be ``None``; all columns must be present.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one taxa axis column is assumed.
            If ``Integral``, then one taxa axis column is assumed.
            If ``None``, then raise error.
        
        taxa_grp_colnames : bool, str, Integral, Sequence of str or None, None, default = True
            Sequence of column names from which to read taxa group axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one taxa group axis column is assumed.
            If ``Integral``, then one taxa group axis column is assumed.
            If ``None``, then raise error.

        trait_colnames : bool, str, Integral, Sequence of str or Integral, None, default = True
            Sequence of column names or indices from which to read trait columns.
            There must be one column name for each trait.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one column name is assumed.
            If ``Integral``, then one column name is assumed.
            If ``None``, then raise error.
        
        value_colname : str, Integral, default = "value"
            Name or index of the value column.
        
        ntaxaaxes : Integral
            Expected number of taxa axes for the matrix.
        
        kwargs : dict
            Additional keyword arguments to be passed for dataframe construction.

        Returns
        -------
        out : pandas.DataFrame
            Output pandas dataframe.
        """
        ###
        ### Type checks and preprocessing
        ###

        # get the number of taxa axes
        check_is_Integral(ntaxaaxes, "ntaxaaxes")

        # get the number of trait axes
        ntraitaxes = 1

        ### process inputs: df
        check_is_pandas_DataFrame(df, "df")

        ### process inputs: taxa_colnames

        # process ``None``
        if taxa_colnames is None:
            taxa_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(taxa_colnames, bool):
            if taxa_colnames:
                taxa_colnames = ["taxa_"+str(i).zfill(len(str(ntaxaaxes))) for i in range(ntaxaaxes)]
            else:
                raise ValueError("must provide taxa column name or index values")
        elif isinstance(taxa_colnames, (str,Integral)):
            taxa_colnames = [taxa_colnames]
        elif isinstance(taxa_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``taxa_colnames`` must be of type ``bool``, ``str``, ``Integral``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(taxa_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(taxa_colnames, "taxa_colnames", ntaxaaxes)

        # check sequence element types
        check_Sequence_all_type(taxa_colnames, "taxa_colnames", (str,Integral))

        # convert taxa_colnames to list so it is mutable
        taxa_colnames = list(taxa_colnames)

        # convert taxa_colnames sequence to all Integral
        for i,elem in enumerate(taxa_colnames):
            if isinstance(elem, str):
                taxa_colnames[i] = df.columns.get_loc(elem)

        # check column indices
        check_pandas_DataFrame_has_column_indices(df, "df", *taxa_colnames)

        ### process inputs: taxa_grp_colnames

        # process ``None``
        if taxa_grp_colnames is None:
            taxa_grp_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(taxa_grp_colnames, bool):
            taxa_grp_colnames = ["taxa_grp_"+str(i).zfill(len(str(ntaxaaxes))) for i in range(ntaxaaxes)] if taxa_grp_colnames else [None for i in range(ntaxaaxes)]
        elif isinstance(taxa_grp_colnames, (str,Integral)):
            taxa_grp_colnames = [taxa_grp_colnames]
        elif isinstance(taxa_grp_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``taxa_grp_colnames`` must be of type ``bool``, ``str``, ``Integral``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(taxa_grp_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(taxa_grp_colnames, "taxa_grp_colnames", ntaxaaxes)

        # check sequence element types
        check_Sequence_all_type(taxa_grp_colnames, "taxa_grp_colnames", (str,Integral,type(None)))

        # convert taxa_grp_colnames to list so it is mutable
        taxa_grp_colnames = list(taxa_grp_colnames)

        # convert taxa_grp_colnames sequence to all Integral
        for i,elem in enumerate(taxa_grp_colnames):
            if isinstance(elem, str):
                taxa_grp_colnames[i] = df.columns.get_loc(elem)

        # create a reduced set that does not contain None (skip None columns)
        reduced = [e for e in taxa_grp_colnames if e is not None]

        # check column indices
        check_pandas_DataFrame_has_column_indices(df, "df", *reduced)

        ### process inputs: trait_colnames

        # process ``None``
        if trait_colnames is None:
            trait_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(trait_colnames, bool):
            if trait_colnames:
                trait_colnames = ["trait_"+str(i).zfill(len(str(ntraitaxes))) for i in range(ntraitaxes)]
            else:
                raise ValueError("must provide trait column name or index values")
        elif isinstance(trait_colnames, (str,Integral)):
            trait_colnames = [trait_colnames]
        elif isinstance(trait_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``trait_colnames`` must be of type ``bool``, ``str``, ``Integral``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(trait_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(trait_colnames, "trait_colnames", ntraitaxes)

        # check sequence element types
        check_Sequence_all_type(trait_colnames, "trait_colnames", (str,type(None)))

        # convert trait_colnames to list so it is mutable
        trait_colnames = list(trait_colnames)

        # convert trait_colnames sequence to all Integral
        for i,elem in enumerate(trait_colnames):
            if isinstance(elem, str):
                trait_colnames[i] = df.columns.get_loc(elem)

        # check column indices
        check_pandas_DataFrame_has_column_indices(df, "df", *trait_colnames)

        ### process inputs: value_colname
        if isinstance(value_colname, str):
            value_colname = df.columns.get_loc(value_colname)
        
        check_is_Integral(value_colname, "value_colname")

        #####
        ##### Processing
        #####

        ###
        ### Process taxa
        ###

        # get taxa columns data
        taxa_data_ls = [df.iloc[:,ix].to_numpy(dtype=object) for ix in taxa_colnames]

        # get unique taxa and corresponding indices
        taxa_vector_ls, taxa_index_ls = tuple(zip(*[numpy.unique(e, return_index=True) for e in taxa_data_ls]))

        # get unique taxa for all columns: union all taxa together
        taxa = reduce(numpy.union1d, taxa_vector_ls)

        # allocate index arrays for each taxa column
        taxaix_ls = [numpy.empty(len(e), dtype=int) for e in taxa_data_ls]

        # calculate taxa indices
        for i,taxon in enumerate(taxa):
            for taxaix,taxa_data in zip(taxaix_ls,taxa_data_ls):
                taxaix[taxa_data == taxon] = i
        
        ###
        ### Process taxa_grp
        ###

        # get taxa_grp columns data
        taxa_grp_data_ls = [None if ix is None else df.iloc[:,ix].to_numpy(dtype=int) for ix in taxa_grp_colnames]

        # get optional taxa group data
        taxa_grp = None
        for taxa_vector,taxa_index,taxa_grp_data in zip(taxa_vector_ls,taxa_index_ls,taxa_grp_data_ls):
            if taxa_grp_data is not None:
                if taxa_grp is None:
                    taxa_grp = numpy.empty(len(taxa), dtype=int)
                for i, taxon in enumerate(taxa):
                    if taxon in taxa_vector:
                        taxa_grp[i] = taxa_grp_data[(taxa_index[taxa_vector == taxon][0])]

        ###
        ### Process trait
        ###

        # get trait column data
        trait_data = df.iloc[:,(trait_colnames[0])].to_numpy(dtype=object)

        # calculate unique 
        trait, traitix = numpy.unique(trait_data, return_inverse=True)

        # combine trait names
        # trait = reduce(numpy.union1d, blah)

        ###
        ### Process values
        ###

        # get value column data
        value_data = df.iloc[:,(value_colname)].to_numpy(dtype=float)

        # get array dimensions
        ntaxa_tup = tuple(len(taxa) for _ in range(ntaxaaxes))
        ntrait_tup = (len(trait),)
        dim_tup = ntaxa_tup + ntrait_tup

        # allocate NaN array for matrix
        mat = numpy.full(dim_tup, numpy.nan, dtype = float)

        # construct index tuple
        ix_tup = tuple(taxaix_ls) + (traitix,)

        # overwirte NaN values with values
        mat[ix_tup] = value_data

        # construct an object
        out = cls(
            mat = mat, 
            taxa = taxa, 
            taxa_grp = taxa_grp, 
            trait = trait, 
        )

        return out

    @classmethod
    def from_csv(
            cls,
            filename: str, 
            taxa_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            taxa_grp_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            trait_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            value_colname: Union[str,Integral] = "value",
            ntaxaaxes: Integral = 2,
            sep: str = ',',
            header: int = 0,
            **kwargs: dict
        ) -> 'DenseSquareBreedingValueMatrix':
        """
        Read a DenseSquareBreedingValueMatrix from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.

        taxa_colnames : bool, str, Integral, Sequence of str or Integral, None, default = True
            Sequence of column names or indices from which to read taxa axis columns.
            There must be one column name or index for each taxa axis.
            If an element in the sequence cannot be ``None``; all columns must be present.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one taxa axis column is assumed.
            If ``Integral``, then one taxa axis column is assumed.
            If ``None``, then raise error.
        
        taxa_grp_colnames : bool, str, Integral, Sequence of str or None, None, default = True
            Sequence of column names from which to read taxa group axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one taxa group axis column is assumed.
            If ``Integral``, then one taxa group axis column is assumed.
            If ``None``, then raise error.

        trait_colnames : bool, str, Integral, Sequence of str or Integral, None, default = True
            Sequence of column names or indices from which to read trait columns.
            There must be one column name for each trait.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one column name is assumed.
            If ``Integral``, then one column name is assumed.
            If ``None``, then raise error.
        
        value_colname : str, Integral, default = "value"
            Name or index of the value column.
        
        ntaxaaxes : Integral
            Expected number of taxa axes for the matrix.

        sep : str
            CSV file separator.
        
        header : int
            Header row index.
        
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : CSVInputOutput
            An object read from a CSV file.
        """
        # read dataframe from file to pandas
        df = pandas.read_csv(
            filepath_or_buffer = filename,
            sep = sep,
            header = header,
            **kwargs
        )

        # convert pandas to matrix
        out = cls.from_pandas(
            df = df,
            taxa_colnames = taxa_colnames,
            taxa_grp_colnames = taxa_grp_colnames,
            trait_colnames = trait_colnames,
            value_colname = value_colname,
            ntaxaaxes = ntaxaaxes,
        )

        return out



def check_is_DenseSquareBreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``DenseSquareBreedingValueMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseSquareBreedingValueMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                DenseSquareBreedingValueMatrix.__name__,
                type(v).__name__
            )
        )
