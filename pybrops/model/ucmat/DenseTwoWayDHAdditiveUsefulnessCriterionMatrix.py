"""
Module implementing dense usefulness criterion matrices.
"""

__all__ = [
    "DenseTwoWayDHAdditiveUsefulnessCriterionMatrix",
    "check_is_DenseTwoWayDHAdditiveUsefulnessCriterionMatrix",
]

import copy
from numbers import Integral, Real
from typing import Optional, Union
import numpy
from scipy.stats import norm

from pybrops.core.error.error_type_numpy import check_is_Real_or_ndarray, check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len, check_ndarray_is_square_along_axes, check_ndarray_ndim
from pybrops.core.mat.DenseSquareTaxaTraitMatrix import DenseSquareTaxaTraitMatrix
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.model.pmgebvmat.DenseTwoWayParentalMeanGenomicEstimatedBreedingValueMatrix import DenseTwoWayParentalMeanGenomicEstimatedBreedingValueMatrix
from pybrops.model.ucmat.UsefulnessCriterionMatrix import UsefulnessCriterionMatrix
from pybrops.model.vmat.DenseTwoWayDHAdditiveGeneticVarianceMatrix import DenseTwoWayDHAdditiveGeneticVarianceMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix


class DenseTwoWayDHAdditiveUsefulnessCriterionMatrix(
        DenseSquareTaxaTraitMatrix,
        UsefulnessCriterionMatrix,
    ):
    """
    docstring for DenseTwoWayDHAdditiveUsefulnessCriterionMatrix.
    """

    ########################## Special Object Methods ##########################

    ##################### Constructor ######################
    def __init__(
            self,
            mat: numpy.ndarray, 
            upper_percentile: Union[numpy.ndarray,Real] = 0.10,
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseTwoWayDHAdditiveUsefulnessCriterionMatrix.
        
        Parameters
        ----------
        mat : numpy.ndarray
            An array of shape ``(n,n,t)`` used to construct the object.

        upper_percentile : numpy.ndarray, Real
            An array of selection upper percentiles of shape ``(t,)``.
            If ``Real``, broadcast to ``(t,)`` array.

        taxa : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.

        taxa_grp : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.

        trait : numpy.ndarray, None
            A numpy.ndarray of shape ``(t,)`` containing trait names.
            If ``None``, do not store any trait name information.

        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # call DenseSquareTaxaTraitMatrix
        super(DenseTwoWayDHAdditiveUsefulnessCriterionMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )
        # must go after constructor since there is a call to self.ntrait
        # process upper_percentile
        if isinstance(upper_percentile, Real):
            upper_percentile = numpy.repeat(upper_percentile, self.ntrait)
        elif isinstance(upper_percentile, numpy.ndarray):
            check_ndarray_ndim(upper_percentile, "upper_percentile", 1)
            check_ndarray_axis_len(upper_percentile, "upper_percentile", 0, self.ntrait)
        else:
            check_is_Real_or_ndarray(upper_percentile, "upper_percentile")
        self._upper_percentile = upper_percentile

    ############## Forward numeric operators ###############
    ### __add__                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### __sub__                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### __mul__                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### __matmul__              inherited from ``DenseSquareTaxaTraitMatrix``
    ### __truediv__             inherited from ``DenseSquareTaxaTraitMatrix``
    ### __floordiv__            inherited from ``DenseSquareTaxaTraitMatrix``
    ### __mod__                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### __divmod__              inherited from ``DenseSquareTaxaTraitMatrix``
    ### __pow__                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### __lshift__              inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rshift__              inherited from ``DenseSquareTaxaTraitMatrix``
    ### __and__                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### __xor__                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### __or__                  inherited from ``DenseSquareTaxaTraitMatrix``

    ############# Backwards numeric operators ##############
    ### __radd__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rsub__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rmul__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rmatmul__             inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rtruediv__            inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rfloordiv__           inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rmod__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rdivmod__             inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rlshift__             inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rrshift__             inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rand__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __rxor__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __ror__                 inherited from ``DenseSquareTaxaTraitMatrix``

    ############# Augmented numeric operators ##############
    ### __iadd__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __isub__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __imul__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __imatmul__             inherited from ``DenseSquareTaxaTraitMatrix``
    ### __itruediv__            inherited from ``DenseSquareTaxaTraitMatrix``
    ### __ifloordiv__           inherited from ``DenseSquareTaxaTraitMatrix``
    ### __imod__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __ipow__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __ilshift__             inherited from ``DenseSquareTaxaTraitMatrix``
    ### __irshift__             inherited from ``DenseSquareTaxaTraitMatrix``
    ### __iand__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __ixor__                inherited from ``DenseSquareTaxaTraitMatrix``
    ### __ior__                 inherited from ``DenseSquareTaxaTraitMatrix``

    ################## Logical operators ###################
    ### __lt__                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### __le__                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### __eq__                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### __ne__                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### __gt__                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### __ge__                  inherited from ``DenseSquareTaxaTraitMatrix``

    ################# Container operators ##################
    ### __len__                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### __getitem__             inherited from ``DenseSquareTaxaTraitMatrix``
    ### __setitem__             inherited from ``DenseSquareTaxaTraitMatrix``
    ### __delitem__             inherited from ``DenseSquareTaxaTraitMatrix``
    ### __iter__                inherited from ``DenseSquareTaxaTraitMatrix``

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseTwoWayDHAdditiveUsefulnessCriterionMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseTwoWayDHAdditiveUsefulnessCriterionMatrix
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
            trait = copy.copy(self.trait),
        )

        # copy taxa metadata
        out.taxa_grp_name = copy.copy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.copy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.copy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.copy(self.taxa_grp_len)

        return out

    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'DenseTwoWayDHAdditiveUsefulnessCriterionMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseTwoWayDHAdditiveUsefulnessCriterionMatrix
        """
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo),
            trait = copy.deepcopy(self.trait, memo)
        )

        # copy taxa metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name, memo)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix, memo)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix, memo)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len, memo)

        return out

    ########### Miscellaneous special functions ############
    ### __repr__                inherited from ``DenseSquareTaxaTraitMatrix``

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    ### mat                     inherited from ``DenseSquareTaxaTraitMatrix``
    @DenseSquareTaxaTraitMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set foo."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim(value, "mat", 3)
        check_ndarray_is_square_along_axes(value, "mat", (0,1))
        self._mat = value

    ############## Matrix Metadata Properties ##############
    ### mat_ndim                inherited from ``DenseSquareTaxaTraitMatrix``
    ### mat_shape               inherited from ``DenseSquareTaxaTraitMatrix``

    ############## Square Metadata Properties ##############
    ### nsquare                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### square_axes             inherited from ``DenseSquareTaxaTraitMatrix``
    ### square_axes_len         inherited from ``DenseSquareTaxaTraitMatrix``

    ########### Square Taxa Metadata Properties ############
    ### nsquare_taxa            inherited from ``DenseSquareTaxaTraitMatrix``
    ### square_taxa_axes        inherited from ``DenseSquareTaxaTraitMatrix``
    ### square_taxa_axes_len    inherited from ``DenseSquareTaxaTraitMatrix``

    ################# Taxa Data Properites #################
    ### taxa                    inherited from ``DenseSquareTaxaTraitMatrix``
    ### taxa_grp                inherited from ``DenseSquareTaxaTraitMatrix``

    ############### Taxa Metadata Properites ###############
    ### ntaxa                   inherited from ``DenseSquareTaxaTraitMatrix``
    ### taxa_axis               inherited from ``DenseSquareTaxaTraitMatrix``
    ### taxa_grp_name           inherited from ``DenseSquareTaxaTraitMatrix``
    ### taxa_grp_stix           inherited from ``DenseSquareTaxaTraitMatrix``
    ### taxa_grp_spix           inherited from ``DenseSquareTaxaTraitMatrix``
    ### taxa_grp_len            inherited from ``DenseSquareTaxaTraitMatrix``

    ###################### Trait data ######################
    ### trait                   inherited from ``DenseSquareTaxaTraitMatrix``

    #################### Trait metadata ####################
    ### ntrait                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### trait_axis              inherited from ``DenseSquareTaxaTraitMatrix``

    ############## Square Metadata Properties ##############
    @DenseSquareTaxaTraitMatrix.square_taxa_axes.getter
    def square_taxa_axes(self) -> tuple:
        """Axis indices for taxa axes that are square."""
        # square taxa axes are everything except the last axis
        return (0,1)

    #################### Trait metadata ####################
    @DenseSquareTaxaTraitMatrix.trait_axis.getter
    def trait_axis(self) -> int:
        """Axis along which traits are stored."""
        return 2

    ################ EPGC metadata property ################
    @property
    def epgc(self) -> tuple:
        """Expected parental genomic contribution to the offspring from each parent."""
        return (0.5,0.5)

    ########### Selection Percentile Properties ############
    @property
    def upper_percentile(self) -> numpy.ndarray:
        """upper_percentile."""
        return self._upper_percentile
    
    @property
    def selection_intensity(self) -> numpy.ndarray:
        """selection_intensity."""
        return norm.pdf(norm.ppf(1.0 - self._upper_percentile)) / self._upper_percentile
    
    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseTwoWayDHAdditiveUsefulnessCriterionMatrix':
        """
        Make a shallow copy of the DenseTwoWayDHAdditiveUsefulnessCriterionMatrix.

        Returns
        -------
        out : DenseTwoWayDHAdditiveUsefulnessCriterionMatrix
            A shallow copy of the original DenseTwoWayDHAdditiveUsefulnessCriterionMatrix.
        """
        return self.__copy__()

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseTwoWayDHAdditiveUsefulnessCriterionMatrix':
        """
        Make a deep copy of the DenseTwoWayDHAdditiveUsefulnessCriterionMatrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseTwoWayDHAdditiveUsefulnessCriterionMatrix
            A deep copy of the original DenseTwoWayDHAdditiveUsefulnessCriterionMatrix.
        """
        return self.__deepcopy__(memo)

    ######### Matrix element copy-on-manipulation ##########
    ### adjoin                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### delete                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### insert                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### select                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### adjoin_taxa             inherited from ``DenseSquareTaxaTraitMatrix``
    ### delete_taxa             inherited from ``DenseSquareTaxaTraitMatrix``
    ### insert_taxa             inherited from ``DenseSquareTaxaTraitMatrix``
    ### select_taxa             inherited from ``DenseSquareTaxaTraitMatrix``
    ### adjoin_trait            inherited from ``DenseSquareTaxaTraitMatrix``
    ### delete_trait            inherited from ``DenseSquareTaxaTraitMatrix``
    ### insert_trait            inherited from ``DenseSquareTaxaTraitMatrix``
    ### select_trait            inherited from ``DenseSquareTaxaTraitMatrix``

    ######### Matrix element in-place-manipulation #########
    ### append                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### remove                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### incorp                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### append_taxa             inherited from ``DenseSquareTaxaTraitMatrix``
    ### remove_taxa             inherited from ``DenseSquareTaxaTraitMatrix``
    ### incorp_taxa             inherited from ``DenseSquareTaxaTraitMatrix``
    ### append_trait            inherited from ``DenseSquareTaxaTraitMatrix``
    ### remove_trait            inherited from ``DenseSquareTaxaTraitMatrix``
    ### incorp_trait            inherited from ``DenseSquareTaxaTraitMatrix``

    ################### Sorting Methods ####################
    ### lexsort                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### reorder                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### sort                    inherited from ``DenseSquareTaxaTraitMatrix``
    ### lexsort_taxa            inherited from ``DenseSquareTaxaTraitMatrix``
    ### reorder_taxa            inherited from ``DenseSquareTaxaTraitMatrix``
    ### sort_taxa               inherited from ``DenseSquareTaxaTraitMatrix``
    ### lexsort_trait           inherited from ``DenseSquareTaxaTraitMatrix``
    ### reorder_trait           inherited from ``DenseSquareTaxaTraitMatrix``
    ### sort_trait              inherited from ``DenseSquareTaxaTraitMatrix``

    ################### Grouping Methods ###################
    ### group                   inherited from ``DenseSquareTaxaTraitMatrix``
    ### ungroup                 inherited from ``DenseSquareTaxaTraitMatrix``
    ### is_grouped              inherited from ``DenseSquareTaxaTraitMatrix``
    ### group_taxa              inherited from ``DenseSquareTaxaTraitMatrix``
    ### ungroup_taxa            inherited from ``DenseSquareTaxaTraitMatrix``
    ### is_grouped_taxa         inherited from ``DenseSquareTaxaTraitMatrix``

    #################### Square Methods ####################
    ### is_square               inherited from ``DenseSquareTaxaTraitMatrix``
    ### is_square_taxa          inherited from ``DenseSquareTaxaTraitMatrix``

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat                  inherited from ``DenseSquareTaxaTraitMatrix``
    ### concat_taxa             inherited from ``DenseSquareTaxaTraitMatrix``
    ### concat_trait            inherited from ``DenseSquareTaxaTraitMatrix``

    ################# Construction Methods #################
    @classmethod
    def from_gmod(
            cls, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nmating: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            upper_percentile: Union[numpy.ndarray,Real],
            **kwargs: dict
        ) -> 'DenseTwoWayDHAdditiveUsefulnessCriterionMatrix':
        """
        Estimate usefulness criterion values from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        nmating : int
            Number of cross patterns to simulate for usefulness criterion
            estimation.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genetic
            covariance.
        nself : int
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.
        gmapfn : GeneticMapFunction
            Genetic map function with which to calculate recombination probabilities.
        upper_percentile : numpy.ndarray, Real
            An array of selection upper percentiles of shape ``(t,)``.
            If ``Real``, broadcast to ``(t,)`` array.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTwoWayDHAdditiveUsefulnessCriterionMatrix
            A matrix of usefulness criterion estimations.
        """
        # check types and inputs
        check_is_GenomicModel(gmod, "gmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")

        # process upper_percentile
        if isinstance(upper_percentile, Real):
            upper_percentile = numpy.repeat(upper_percentile, gmod.ntrait)
        elif isinstance(upper_percentile, numpy.ndarray):
            check_ndarray_ndim(upper_percentile, "upper_percentile", 1)
            check_ndarray_axis_len(upper_percentile, "upper_percentile", 0, gmod.ntrait)
        else:
            check_is_Real_or_ndarray(upper_percentile, "upper_percentile")

        # create dense pmGEBV matrix
        # (n,n,t)
        meanmat = DenseTwoWayParentalMeanGenomicEstimatedBreedingValueMatrix.from_gmod(
            gmod = gmod, 
            gmat = pgmat
        )

        # create variance matrix
        # (n,n,t)
        varmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_gmod(
            gmod = gmod,
            pgmat = pgmat,
            nmating = nmating,
            nprogeny = nprogeny,
            nself = nself,
            gmapfn = gmapfn,
        )

        # get cross means
        # (n,n,t)
        p_mean = meanmat.unscale()

        # get cross variances
        # (n,n,t)
        p_var = varmat.mat

        # calculate cross standard deviations
        # (n,n,t)
        p_std = numpy.sqrt(p_var)

        # calculate selection intensity
        # (1,1,t,)
        p_si = norm.pdf(norm.ppf(1.0 - upper_percentile[None,None,:])) / upper_percentile[None,None,:]

        # calculate usefulness criterion
        # (n,n,t) + (1,1,t) * (n,n,t)
        p_uc = p_mean + p_si * p_std

        # construct object
        out = cls(
            mat = p_uc,
            upper_percentile = upper_percentile,
            taxa = meanmat.taxa,
            taxa_grp = meanmat.taxa_grp,
            trait = meanmat.trait,
        )

        # copy metadata if any
        out.taxa_grp_len  = copy.copy(meanmat.taxa_grp_len)
        out.taxa_grp_name = copy.copy(meanmat.taxa_grp_name)
        out.taxa_grp_stix = copy.copy(meanmat.taxa_grp_stix)
        out.taxa_grp_spix = copy.copy(meanmat.taxa_grp_spix)

        return out

    @classmethod
    def from_algmod(
            cls, 
            algmod: AdditiveLinearGenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nmating: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            upper_percentile: Union[numpy.ndarray,Real],
            mem: Union[Integral,None] = 1024,
            **kwargs: dict
        ) -> 'DenseTwoWayDHAdditiveUsefulnessCriterionMatrix':
        """
        Estimate usefulness criterion values from an AdditiveLinearGenomicModel.

        Parameters
        ----------
        algmod : AdditiveLinearGenomicModel
            AdditiveLinearGenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        nmating : int
            Number of cross patterns to simulate for usefulness criterion
            estimation.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genetic
            covariance.
        nself : int
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.
        gmapfn : GeneticMapFunction
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        upper_percentile : numpy.ndarray, Real
            An array of selection upper percentiles of shape ``(t,)``.
            If ``Real``, broadcast to ``(t,)`` array.
        mem : int
            Memory chunk size to use during matrix operations.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTwoWayDHAdditiveUsefulnessCriterionMatrix
            A matrix of additive usefulness criterion estimations.
        """
        # check types and inputs
        check_is_GenomicModel(algmod, "algmod")
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")

        # process upper_percentile
        if isinstance(upper_percentile, Real):
            upper_percentile = numpy.repeat(upper_percentile, algmod.ntrait)
        elif isinstance(upper_percentile, numpy.ndarray):
            check_ndarray_ndim(upper_percentile, "upper_percentile", 1)
            check_ndarray_axis_len(upper_percentile, "upper_percentile", 0, algmod.ntrait)
        else:
            check_is_Real_or_ndarray(upper_percentile, "upper_percentile")

        # create dense pmGEBV matrix
        # (n,n,t)
        meanmat = DenseTwoWayParentalMeanGenomicEstimatedBreedingValueMatrix.from_gmod(
            gmod = algmod, 
            gmat = pgmat
        )

        # create variance matrix
        varmat = DenseTwoWayDHAdditiveGeneticVarianceMatrix.from_algmod(
            algmod = algmod,
            pgmat = pgmat,
            nmating = nmating,
            nprogeny = nprogeny,
            nself = nself,
            gmapfn = gmapfn,
            mem = mem,
        )

        # get cross means
        # (n,n,t)
        p_mean = meanmat.unscale()

        # get cross variances
        # (n,n,t)
        p_var = varmat.mat

        # calculate cross standard deviations
        # (n,n,t)
        p_std = numpy.sqrt(p_var)

        # calculate selection intensity
        # (1,1,t,)
        p_si = norm.pdf(norm.ppf(1.0 - upper_percentile[None,None,:])) / upper_percentile[None,None,:]

        # calculate usefulness criterion
        # (n,n,t) + (1,1,t) * (n,n,t)
        p_uc = p_mean + p_si * p_std

        # construct object
        out = cls(
            mat = p_uc,
            upper_percentile = upper_percentile,
            taxa = meanmat.taxa,
            taxa_grp = meanmat.taxa_grp,
            trait = meanmat.trait,
        )

        # copy metadata if any
        out.taxa_grp_len  = copy.copy(meanmat.taxa_grp_len)
        out.taxa_grp_name = copy.copy(meanmat.taxa_grp_name)
        out.taxa_grp_stix = copy.copy(meanmat.taxa_grp_stix)
        out.taxa_grp_spix = copy.copy(meanmat.taxa_grp_spix)

        return out



################################## Utilities ###################################
def check_is_DenseTwoWayDHAdditiveUsefulnessCriterionMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``DenseTwoWayDHAdditiveUsefulnessCriterionMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseTwoWayDHAdditiveUsefulnessCriterionMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                DenseTwoWayDHAdditiveUsefulnessCriterionMatrix.__name__,
                type(v).__name__
            )
        )
