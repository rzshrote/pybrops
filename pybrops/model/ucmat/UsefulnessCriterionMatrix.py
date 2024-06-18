"""
Module defining interfaces and error checking routines for usefulness criterion value matrices.
"""

__all__ = [
    "UsefulnessCriterionMatrix",
    "check_is_UsefulnessCriterionMatrix",
]

from abc import ABCMeta, abstractmethod

import numpy
from pybrops.core.mat.SquareTaxaTraitMatrix import SquareTaxaTraitMatrix
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class UsefulnessCriterionMatrix(
        SquareTaxaTraitMatrix,
        metaclass = ABCMeta,
    ):
    """
    Abstract class for usefulness criterion matrix representation.
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

    ################ EPGC metadata property ################
    @property
    @abstractmethod
    def epgc(self) -> tuple:
        """Expected parental genomic contribution to the offspring from each parent."""
        raise NotImplementedError("property is abstract")

    ########### Selection Percentile Properties ############
    @property
    @abstractmethod
    def upper_percentile(self) -> numpy.ndarray:
        """Upper percentile of selection."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def selection_intensity(self) -> numpy.ndarray:
        """Corresponding to the upper percentile of selection."""
        raise NotImplementedError("property is abstract")

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

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat                  inherited from ``SquareTaxaTraitMatrix``
    ### concat_taxa             inherited from ``SquareTaxaTraitMatrix``
    ### concat_trait            inherited from ``SquareTaxaTraitMatrix``

    ################# Construction Methods #################
    @classmethod
    @abstractmethod
    def from_gmod(
            cls, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nmating: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            **kwargs: dict
        ) -> 'UsefulnessCriterionMatrix':
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : UsefulnessCriterionMatrix
            A matrix of usefulness criterion estimations.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_UsefulnessCriterionMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``UsefulnessCriterionMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, UsefulnessCriterionMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                UsefulnessCriterionMatrix.__name__,
                type(v).__name__
            )
        )
