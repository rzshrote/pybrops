"""
Module implementing representation of dense two-way progeny mean EBV matrices.
"""

__all__ = [
    "DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix",
    "check_is_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix",
]

from numbers import Integral, Real
from typing import Optional, Union
import numpy
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_all_gteq, check_ndarray_axis_len, check_ndarray_is_square, check_ndarray_ndim
from pybrops.core.error.error_value_python import check_is_gteq
from pybrops.core.mat.DenseSquareTaxaTraitMatrix import DenseSquareTaxaTraitMatrix
from pybrops.model.pmebvmat.ProgenyMeanEstimatedBreedingValueMatrix import ProgenyMeanEstimatedBreedingValueMatrix


class DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix(
        DenseSquareTaxaTraitMatrix,
        ProgenyMeanEstimatedBreedingValueMatrix,
    ):
    """
    A concrete class for dense progeny mean EBV matrix representation.

    The purpose of this concrete class is to implement functionality for:
        1) Progeny mean EBV estimation from breeding value matrices.
        2) I/O for two-way progeny mean EBV matrices.
    """

    ########################## Special Object Methods ##########################
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
        Constructor for DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix.
        
        Parameters
        ----------
        mat : numpy.ndarray
            An array of progeny mean estimated breeding values of shape 
            ``(ntaxa,ntaxa,ntrait)``.

            Where:

            - ``ntaxa`` is the number of taxa.
            - ``ntrait`` is the number of traits.

            It is the responsibility of the user to ensure that the means and 
            standard deviations of this array along the ``taxa`` axis are ``0`` and
            ``1``, respectively, if the breeding values are with respect to the
            individuals in the breeding value matrix.

        location : numpy.ndarray, Real
            If ``numpy.ndarray``, an array of shape ``(ntrait,)`` containing 
            breeding value locations.
            If ``Real``, create a ``numpy.ndarray`` of shape ``(ntrait,)`` 
            filled with the provided value.
        
        scale : numpy.ndarray, Real
            If ``numpy.ndarray``, an array of shape ``(ntrait,)`` containing 
            breeding value scales.
            If ``Real``, create a ``numpy.ndarray`` of shape ``(ntrait,)`` 
            filled with the provided value.
        
        taxa : numpy.ndarray, None
            If ``numpy.ndarray``, an array of shape ``(ntaxa,)`` containing 
            taxa names.
            If ``None``, do not store any taxa name information.
        
        taxa_grp : numpy.ndarray, None
            If ``numpy.ndarray``, an array of shape ``(ntaxa,)`` containing 
            taxa groupings.
            If ``None``, do not store any taxa group information.
        
        trait : numpy.ndarray, None
            If ``numpy.ndarray``, an array of shape ``(ntrait,)`` containing 
            trait names.
            If ``None``, do not store any trait name information.
        
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # call DenseSquareTaxaTraitMatrix constructor
        super(DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )
        # set location and scale parameters
        self.location = location
        self.scale = scale


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

    ##################### Matrix Data ######################
    @DenseSquareTaxaTraitMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set pointer to raw numpy.ndarray object."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim(value, "mat", 3) # (ntaxa,ntaxa,ntrait)

    ################# Breeding Value Data ##################
    @property
    def location(self) -> numpy.ndarray:
        """Mean of the phenotype values used to calculate breeding values."""
        return self._location
    @location.setter
    def location(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set the mean of the phenotype values used to calculate breeding values"""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "location", 1)
            check_ndarray_axis_len(value, "location", 0, self.ntrait)
        elif isinstance(value, Real):
            value = numpy.repeat(value, self.ntrait)
        else:
            raise TypeError("variable 'location' must be of type 'numpy.ndarray' or 'Real'")
        self._location = value
    
    @property
    def scale(self) -> numpy.ndarray:
        """Standard deviation of the phenotype values used to calculate breeding values."""
        return self._scale
    @scale.setter
    def scale(self, value: Union[numpy.ndarray,Real]) -> None:
        """Set the standard deviation of the phenotype values used to calculate breeding values"""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "scale", 1)
            check_ndarray_axis_len(value, "scale", 0, self.ntrait)
            check_ndarray_all_gteq(value, "scale", 0)
        elif isinstance(value, Real):
            check_is_gteq(value, "scale", 0)
            value = numpy.repeat(value, self.ntrait)
        else:
            raise TypeError("variable 'scale' must be of type 'numpy.ndarray' or 'Real'")
        self._scale = value

    ############## Square Metadata Properties ##############
    @DenseSquareTaxaTraitMatrix.square_taxa_axes.getter
    def square_taxa_axes(self) -> tuple:
        """Axis indices for taxa axes that are square."""
        return (0,1) # (female, male)

    @DenseSquareTaxaTraitMatrix.trait_axis.getter
    def trait_axis(self) -> Integral:
        """Axis along which traits are stored."""
        return 2

    ######## Expected parental genome contributions ########
    @property
    def epgc(self) -> tuple:
        """Expected parental genome contribution to the offspring from each parent."""
        return (0.5, 0.5) # (female, male)    

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

    ############################## Class Methods ###############################

    #################### Import Methods ####################

    @classmethod
    def from_numpy(
            cls, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'ProgenyMeanEstimatedBreedingValueMatrix':
        """
        Construct a ProgenyMeanEstimatedBreedingValueMatrix from a numpy.ndarray.
        Calculates mean-centering and scaling to unit variance.

        Parameters
        ----------
        mat : numpy.ndarray
            An array of shape ``(ntaxa,ntaxa,ntrait)`` for which to mean-center and scale.

            Where:

            - ``ntaxa`` is the number of taxa.
            - ``ntrait`` is the number of traits.

        taxa : numpy.ndarray
            An array of taxa names of shape ``(ntaxa,)``.

        taxa_grp : numpy.ndarray
            An array of taxa groups of shape ``(ntaxa)``.

        trait : numpy.ndarray
            An array of trait names of shape ``(ntrait,)``.

        Returns
        -------
        out : ProgenyMeanEstimatedBreedingValueMatrix
            Output progeny mean estimated breeding value matrix.
        """
        # check type inputs
        check_is_ndarray(mat, "mat")
        check_ndarray_ndim(mat, "mat", 3)

        # calculate location parameters
        # (n,n,t) -> (t,)
        location = numpy.nanmean(mat, axis = (0,1))
        
        # calculate location scale parameters
        # (n,n,t) -> (t,)
        scale = numpy.nanstd(mat, axis = (0,1))

        # if scale == 0.0, then set to 1.0 (do not scale to avoid division by zero)
        scale[scale == 0.0] = 1.0

        # mean center and scale values
        # scalar / (1,t) -> (1,t)
        # (1,t) * ( (n,t) - (1,t) ) -> (n,t)
        # multiply since multiplication is faster than division for floating points
        mat = (1.0 / scale[None,None,:]) * (mat - location[None,None,:]) 

        # construct output
        out = cls(
            mat = mat,
            location = location,
            scale = scale,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

        return out

    ################# Construction Methods #################

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix.__name__,
                type(v).__name__
            )
        )
