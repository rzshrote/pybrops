"""
Module implementing dense two-way parental mean GEBV matrices and associated error checking routines.
"""

__all__ = [
    "DenseTwoWayParentalMeanEstimatedBreedingValueMatrix",
    "check_is_DenseTwoWayParentalMeanEstimatedBreedingValueMatrix",
]

import copy
from itertools import chain
from numbers import Integral
from numbers import Real
from typing import Optional
from typing import Union
import numpy
from pybrops.core.error.error_type_numpy import check_is_Real_or_ndarray
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_bool
from pybrops.core.error.error_value_numpy import check_ndarray_all_gteq
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len
from pybrops.core.error.error_value_numpy import check_ndarray_is_square_along_axes
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_python import check_is_gteq
from pybrops.core.mat.DenseScaledSquareTaxaTraitMatrix import DenseScaledSquareTaxaTraitMatrix
from pybrops.model.pmebvmat.ParentalMeanEstimatedBreedingValueMatrix import ParentalMeanEstimatedBreedingValueMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix, check_is_BreedingValueMatrix

class DenseTwoWayParentalMeanEstimatedBreedingValueMatrix(
        DenseScaledSquareTaxaTraitMatrix,
        ParentalMeanEstimatedBreedingValueMatrix,
    ):
    """
    A concrete class for dense parental mean GEBV matrix representation.
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
        Constructor for DenseTwoWayParentalMeanEstimatedBreedingValueMatrix.
        
        Parameters
        ----------
        mat : numpy.ndarray
            An array of shape ``(n,n,t)`` containing parental mean GEBVs.

            Where:

            - ``n`` is the number of taxa.
            - ``t`` is the number of traits.

        location : numpy.ndarray, Real
            An array of shape ``(t,)`` containing locations for each trait.
            If ``Real``, then the provided location is used for each trait.

        scale : numpy.ndarray, Real
            An array of shape ``(t,)`` containing scales for each trait.
            If ``Real``, then the provided scale is used for each trait.

        taxa : numpy.ndarray, None
            An array of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.

        taxa_grp : numpy.ndarray, None
            An array of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.

        trait : numpy.ndarray, None
            An array of shape ``(t,)`` containing trait names.
            If ``None``, do not store any trait name information.

        kwargs : dict
            Additional keyword arguments.
        """
        # call DenseScaledSquareTaxaTraitMatrix constructor
        super(DenseTwoWayParentalMeanEstimatedBreedingValueMatrix, self).__init__(
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
        ) -> 'DenseTwoWayParentalMeanEstimatedBreedingValueMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseTwoWayParentalMeanEstimatedBreedingValueMatrix
            A copy of the DenseTwoWayParentalMeanEstimatedBreedingValueMatrix.
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
        ) -> 'DenseTwoWayParentalMeanEstimatedBreedingValueMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict, None, default = None
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseTwoWayParentalMeanEstimatedBreedingValueMatrix
            A deep copy of the DenseTwoWayParentalMeanEstimatedBreedingValueMatrix.
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
    def __repr__(
            self
        ) -> str:
        """
        Return repr(self).
        
        Returns
        -------
        out : str
            A representation of the object.
        """
        return "<{0} of shape (nfemale = {1}, nmale = {2}, ntrait = {3}) at {4}>".format(
            type(self).__name__,
            self.nfemale,
            self.nmale,
            self.ntrait,
            hex(id(self)),
        )

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    @DenseScaledSquareTaxaTraitMatrix.mat.setter
    def mat(self, value: object) -> None:
        """Set pointer to raw numpy.ndarray object."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim(value, "mat", 3) # (ntaxa,ntaxa,ntrait)
        check_ndarray_is_square_along_axes(value, "mat", self.square_taxa_axes)
        self._mat = value

    ############## Matrix Metadata Properties ##############
    ### mat_ndim                inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### mat_shape               inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ############## Square Metadata Properties ##############
    ### nsquare                 inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### square_axes             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### square_axes_len         inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ########### Square Taxa Metadata Properties ############
    ### nsquare_taxa            inherited from ``DenseScaledSquareTaxaTraitMatrix``
    
    @property
    def square_taxa_axes(self) -> tuple:
        """square_taxa_axes."""
        return (0,1)
    
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
    
    @property
    def trait_axis(self) -> int:
        """Axis along which traits are stored."""
        return 2

    ############## Scale metadata properties ###############
    @DenseScaledSquareTaxaTraitMatrix.location.setter
    def location(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set the location of the matrix data."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "location", 1)
            check_ndarray_axis_len(value, "location", 0, self.ntrait)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.ntrait)
        else:
            check_is_Real_or_ndarray(value, "location")
        self._location = value
    
    @DenseScaledSquareTaxaTraitMatrix.scale.setter
    def scale(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set the scale of the matrix data."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "scale", 1)
            check_ndarray_axis_len(value, "scale", 0, self.ntrait)
            check_ndarray_all_gteq(value, "scale", 0)
        elif isinstance(value, Real):
            check_is_gteq(value, "scale", 0)
            value = numpy.repeat(value, self.ntrait)
        else:
            check_is_Real_or_ndarray(value, "scale")
        self._scale = value

    ################ EPGC metadata property ################
    @property
    def epgc(self) -> tuple:
        """Expected parental genomic contribution."""
        return (0.5, 0.5) # female, male

    ################# Parental dimensions ##################
    @property
    def nfemale(self) -> Integral:
        """Number of female parents."""
        return self._mat.shape[self.female_axis]
    
    @property
    def female_axis(self) -> Integral:
        """Axis along which female parents are stored."""
        return 0
    
    @property
    def nmale(self) -> Integral:
        """Number of male parents."""
        return self._mat.shape[self.male_axis]
    
    @property
    def male_axis(self) -> Integral:
        """Axis along which male parents are stored."""
        return 1

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseTwoWayParentalMeanEstimatedBreedingValueMatrix':
        """
        Make a shallow copy of the DenseTwoWayParentalMeanEstimatedBreedingValueMatrix.

        Returns
        -------
        out : DenseTwoWayParentalMeanEstimatedBreedingValueMatrix
            A shallow copy of the original DenseTwoWayParentalMeanEstimatedBreedingValueMatrix.
        """
        return self.__copy__()

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseTwoWayParentalMeanEstimatedBreedingValueMatrix':
        """
        Make a deep copy of the DenseTwoWayParentalMeanEstimatedBreedingValueMatrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseTwoWayParentalMeanEstimatedBreedingValueMatrix
            A deep copy of the original DenseTwoWayParentalMeanEstimatedBreedingValueMatrix.
        """
        return self.__deepcopy__(memo)

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
    def transform(self, mat: numpy.ndarray, copy: bool = False) -> numpy.ndarray:
        """
        Transform a provided matrix using location and scale parameters within 
        the matrix.

        Parameters
        ----------
        mat : numpy.ndarray
            An array of shape ``(?,?,t)`` to be centered and scaled.

            Where

            - ``?`` can be any dimension length.
            - ``t`` is the number of traits. Must be equivalent to the number 
              of traits represented by the matrix.

        copy : bool, default = False
            Whether to copy the input ``mat`` or not.

        Returns
        -------
        out : numpy.ndarray
            A transformed array.
        """
        # type checks
        check_is_ndarray(mat, "mat")
        check_ndarray_ndim(mat, "mat", self.mat_ndim)
        check_ndarray_axis_len(mat, "mat", self.trait_axis, self.ntrait)
        check_is_bool(copy, "copy")

        # get dimensions we'll need
        trait_axis = self.trait_axis
        naxes = self.mat_ndim

        # create tuple for reshaping axes for broadcasting
        # (None,...,slice(None),...,None)
        broadcast_shape = tuple(None if i != trait_axis else slice(None) for i in range(naxes))

        # copy matrix if needed
        out = mat.copy() if copy else mat

        # transform values
        # (?,...,t,...,?) -= (1,...,t,...,1)
        # (?,...,t,...,?) *= (1,...,t,...,1)
        out -= self.location[broadcast_shape]
        out *= (1.0 / self.scale[broadcast_shape])

        return out

    def untransform(self, mat: numpy.ndarray, copy: bool = False) -> numpy.ndarray:
        """
        Untransform a provided matrix using location and scale parameters
        within the DenseScaledMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            An array of shape ``(?,?,t)`` to be unscaled and decentered.

            Where

            - ``?`` can be any dimension length.
            - ``t`` is the number of traits. Must be equivalent to the number 
              of traits represented by the matrix.

        copy : bool, default = False
            Whether to copy the input ``mat`` or not.

        Returns
        -------
        out : numpy.ndarray
            A transformed array.
        """
        # type checks
        check_is_ndarray(mat, "mat")
        check_ndarray_ndim(mat, "mat", self.mat_ndim)
        check_ndarray_axis_len(mat, "mat", self.trait_axis, self.ntrait)
        check_is_bool(copy, "copy")

        # get dimensions we'll need
        trait_axis = self.trait_axis
        naxes = self.mat_ndim

        # create tuple for reshaping axes for broadcasting
        # (None,...,slice(None),...,None)
        broadcast_shape = tuple(None if i != trait_axis else slice(None) for i in range(naxes))

        # copy matrix if needed
        out = mat.copy() if copy else mat

        # transform values
        # (?,...,t,...,?) *= (1,...,t,...,1)
        # (?,...,t,...,?) += (1,...,t,...,1)
        out *= self.scale[broadcast_shape]
        out += self.location[broadcast_shape]

        return out
    
    def rescale(self, inplace: bool = True) -> numpy.ndarray:
        """
        Transform matrix values within the matrix to have centered and scaled 
        values.

        Parameters
        ----------
        inplace : bool, default = True
            Whether to modify matrix values in-place.
            If ``True``, then values in the matrix are scaled internally.
            If ``False``, then values in the matrix are not scaled internally 
            and instead a scaled copy of the values is returned.
        
        Returns
        -------
        out : numpy.ndarray
            A matrix of centered and scaled values. If ``inplace == True``, then
            a pointer to the raw matrix is returned.
        """
        # type checks
        check_is_bool(inplace, "inplace")

        # get dimensions we'll need
        trait_axis = self.trait_axis
        naxes = self.mat_ndim

        # create tuple for reshaping axes for broadcasting
        # (None,...,slice(None),...,None)
        broadcast_shape = tuple(None if i != trait_axis else slice(None) for i in range(naxes))

        # get matrix on which to work
        out = self.mat if inplace else self.mat.copy()

        # unscale values
        out *= self.scale[broadcast_shape]
        out += self.location[broadcast_shape]

        # get the axes along which to calculate the location and scale (exclude taxa axis)
        # (0,...,t-1,t+1,...,n)
        axes = tuple(chain(range(trait_axis),range(trait_axis+1,naxes)))

        # calculate new location parameters
        # (?,...,t,...,?) -> (t,)
        new_location = numpy.nanmean(out, axis = axes)

        # calculate new scale parameters
        # (?,...,t,...,?) -> (t,)
        new_scale = numpy.nanstd(out, axis = axes)
        new_scale[new_scale == 0.0] = 1.0   # if scale == 0.0, set to 1.0 (do not scale)

        # rescale values
        # (?,...,t,...,?) -= (1,...,t,...,1)
        # (?,...,t,...,?) *= (1,...,t,...,1)
        out -= new_location[broadcast_shape]
        out *= (1.0 / new_scale[broadcast_shape])

        # if in-place modification, update location, scale
        if inplace:
            self.location = new_location    # set location
            self.scale = new_scale          # set scale

        return out

    def unscale(self, inplace: bool = True) -> numpy.ndarray:
        """
        Transform matrix values within the matrix back to their unscaled 
        and decentered values.

        Parameters
        ----------
        inplace : bool, default = True
            Whether to modify matrix values in-place.
        
        Returns
        -------
        out : numpy.ndarray
            A matrix of unscaled and decentered values. If ``inplace == True``,
            then a pointer to the raw matrix is returned.
        """
        # type checks
        check_is_bool(inplace, "inplace")

        # get dimensions we'll need
        trait_axis = self.trait_axis
        naxes = self.mat_ndim

        # create tuple for reshaping axes for broadcasting
        # (None,...,slice(None),...,None)
        broadcast_shape = tuple(None if i != trait_axis else slice(None) for i in range(naxes))

        # get matrix on which to work
        out = self.mat if inplace else self.mat.copy()

        # unscale and uncenter values
        # (?,...,t,...,?) * (1,...,t,...,1) -> (?,...,t,...,?)
        # (?,...,t,...,?) + (1,...,t,...,1) -> (?,...,t,...,?)
        out *= self.scale[broadcast_shape]
        out += self.location[broadcast_shape]

        # if inplace, set scale and location to 1.0 and 0.0, respectively
        if inplace:
            self.scale[:] = 1.0
            self.location[:] = 0.0

        return out

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat                  inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### concat_taxa             inherited from ``DenseScaledSquareTaxaTraitMatrix``
    ### concat_trait            inherited from ``DenseScaledSquareTaxaTraitMatrix``

    ################# Construction Methods #################
    @classmethod
    def from_bvmat(
            cls,
            bvmat: BreedingValueMatrix,
            **kwargs: dict
        ) -> 'DenseTwoWayParentalMeanEstimatedBreedingValueMatrix':
        """
        Calculate parental mean GEBVs for each pairwise 2-way cross between individuals.

        Parameters
        ----------
        bvmat : BreedingValueMatrix
            Estimated breeding value matrix with which to calculate parental means.

        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTwoWayParentalMeanEstimatedBreedingValueMatrix
            A matrix of parental mean GEBVs for a two way cross.
        """
        # type checks
        check_is_BreedingValueMatrix(bvmat, "bvmat")

        # get parental matrix
        # (n,t)
        pmat = bvmat.mat

        # reshape, add matrices together, and divide by 2 to get parental mean matrix
        # scalar * ((n,1,t) + (1,n,t)) = (n,n,t)
        pmmat = 0.5 * (pmat[:,None,:] + pmat[None,:,:])

        # copy metadata from GEBV matrix and construct matrix
        out = cls(
            mat = pmmat,
            location = bvmat.location,
            scale = bvmat.scale,
            taxa = bvmat.taxa,
            taxa_grp = bvmat.taxa_grp,
            trait = bvmat.trait,
        )

        # copy metadata
        out.taxa_grp_name = bvmat.taxa_grp_name
        out.taxa_grp_stix = bvmat.taxa_grp_stix
        out.taxa_grp_spix = bvmat.taxa_grp_spix
        out.taxa_grp_len  = bvmat.taxa_grp_len

        return out

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_DenseTwoWayParentalMeanEstimatedBreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``DenseTwoWayParentalMeanEstimatedBreedingValueMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseTwoWayParentalMeanEstimatedBreedingValueMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                DenseTwoWayParentalMeanEstimatedBreedingValueMatrix.__name__,
                type(v).__name__
            )
        )
