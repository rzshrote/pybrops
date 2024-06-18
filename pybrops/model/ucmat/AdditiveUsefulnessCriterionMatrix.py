"""
Module defining interfaces and error checking routines for additive usefulness criterion value matrices.
"""

__all__ = [
    "AdditiveUsefulnessCriterionMatrix",
    "check_is_AdditiveUsefulnessCriterionMatrix",
]

from abc import ABCMeta, abstractmethod

from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.ucmat.UsefulnessCriterionMatrix import UsefulnessCriterionMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class AdditiveUsefulnessCriterionMatrix(
        UsefulnessCriterionMatrix,
        metaclass = ABCMeta,
    ):
    """
    Abstract class for additive usefulness criterion matrix representation.
    """
    ########################## Special Object Methods ##########################

    ############## Forward numeric operators ###############
    ### __add__                 inherited from ``UsefulnessCriterionMatrix``
    ### __sub__                 inherited from ``UsefulnessCriterionMatrix``
    ### __mul__                 inherited from ``UsefulnessCriterionMatrix``
    ### __matmul__              inherited from ``UsefulnessCriterionMatrix``
    ### __truediv__             inherited from ``UsefulnessCriterionMatrix``
    ### __floordiv__            inherited from ``UsefulnessCriterionMatrix``
    ### __mod__                 inherited from ``UsefulnessCriterionMatrix``
    ### __divmod__              inherited from ``UsefulnessCriterionMatrix``
    ### __pow__                 inherited from ``UsefulnessCriterionMatrix``
    ### __lshift__              inherited from ``UsefulnessCriterionMatrix``
    ### __rshift__              inherited from ``UsefulnessCriterionMatrix``
    ### __and__                 inherited from ``UsefulnessCriterionMatrix``
    ### __xor__                 inherited from ``UsefulnessCriterionMatrix``
    ### __or__                  inherited from ``UsefulnessCriterionMatrix``

    ############# Backwards numeric operators ##############
    ### __radd__                inherited from ``UsefulnessCriterionMatrix``
    ### __rsub__                inherited from ``UsefulnessCriterionMatrix``
    ### __rmul__                inherited from ``UsefulnessCriterionMatrix``
    ### __rmatmul__             inherited from ``UsefulnessCriterionMatrix``
    ### __rtruediv__            inherited from ``UsefulnessCriterionMatrix``
    ### __rfloordiv__           inherited from ``UsefulnessCriterionMatrix``
    ### __rmod__                inherited from ``UsefulnessCriterionMatrix``
    ### __rdivmod__             inherited from ``UsefulnessCriterionMatrix``
    ### __rlshift__             inherited from ``UsefulnessCriterionMatrix``
    ### __rrshift__             inherited from ``UsefulnessCriterionMatrix``
    ### __rand__                inherited from ``UsefulnessCriterionMatrix``
    ### __rxor__                inherited from ``UsefulnessCriterionMatrix``
    ### __ror__                 inherited from ``UsefulnessCriterionMatrix``

    ############# Augmented numeric operators ##############
    ### __iadd__                inherited from ``UsefulnessCriterionMatrix``
    ### __isub__                inherited from ``UsefulnessCriterionMatrix``
    ### __imul__                inherited from ``UsefulnessCriterionMatrix``
    ### __imatmul__             inherited from ``UsefulnessCriterionMatrix``
    ### __itruediv__            inherited from ``UsefulnessCriterionMatrix``
    ### __ifloordiv__           inherited from ``UsefulnessCriterionMatrix``
    ### __imod__                inherited from ``UsefulnessCriterionMatrix``
    ### __ipow__                inherited from ``UsefulnessCriterionMatrix``
    ### __ilshift__             inherited from ``UsefulnessCriterionMatrix``
    ### __irshift__             inherited from ``UsefulnessCriterionMatrix``
    ### __iand__                inherited from ``UsefulnessCriterionMatrix``
    ### __ixor__                inherited from ``UsefulnessCriterionMatrix``
    ### __ior__                 inherited from ``UsefulnessCriterionMatrix``

    ################## Logical operators ###################
    ### __lt__                  inherited from ``UsefulnessCriterionMatrix``
    ### __le__                  inherited from ``UsefulnessCriterionMatrix``
    ### __eq__                  inherited from ``UsefulnessCriterionMatrix``
    ### __ne__                  inherited from ``UsefulnessCriterionMatrix``
    ### __gt__                  inherited from ``UsefulnessCriterionMatrix``
    ### __ge__                  inherited from ``UsefulnessCriterionMatrix``

    ################# Container operators ##################
    ### __len__                 inherited from ``UsefulnessCriterionMatrix``
    ### __getitem__             inherited from ``UsefulnessCriterionMatrix``
    ### __setitem__             inherited from ``UsefulnessCriterionMatrix``
    ### __delitem__             inherited from ``UsefulnessCriterionMatrix``
    ### __iter__                inherited from ``UsefulnessCriterionMatrix``

    #################### Matrix copying ####################
    ### __copy__                inherited from ``UsefulnessCriterionMatrix``
    ### __deepcopy__            inherited from ``UsefulnessCriterionMatrix``

    ########### Miscellaneous special functions ############
    ### __repr__                inherited from ``UsefulnessCriterionMatrix``

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    ### mat                     inherited from ``UsefulnessCriterionMatrix``

    ############## Matrix Metadata Properties ##############
    ### mat_ndim                inherited from ``UsefulnessCriterionMatrix``
    ### mat_shape               inherited from ``UsefulnessCriterionMatrix``

    ############## Square Metadata Properties ##############
    ### nsquare                 inherited from ``UsefulnessCriterionMatrix``
    ### square_axes             inherited from ``UsefulnessCriterionMatrix``
    ### square_axes_len         inherited from ``UsefulnessCriterionMatrix``

    ########### Square Taxa Metadata Properties ############
    ### nsquare_taxa            inherited from ``UsefulnessCriterionMatrix``
    ### square_taxa_axes        inherited from ``UsefulnessCriterionMatrix``
    ### square_taxa_axes_len    inherited from ``UsefulnessCriterionMatrix``

    ################# Taxa Data Properites #################
    ### taxa                    inherited from ``UsefulnessCriterionMatrix``
    ### taxa_grp                inherited from ``UsefulnessCriterionMatrix``

    ############### Taxa Metadata Properites ###############
    ### ntaxa                   inherited from ``UsefulnessCriterionMatrix``
    ### taxa_axis               inherited from ``UsefulnessCriterionMatrix``
    ### taxa_grp_name           inherited from ``UsefulnessCriterionMatrix``
    ### taxa_grp_stix           inherited from ``UsefulnessCriterionMatrix``
    ### taxa_grp_spix           inherited from ``UsefulnessCriterionMatrix``
    ### taxa_grp_len            inherited from ``UsefulnessCriterionMatrix``

    ###################### Trait data ######################
    ### trait                   inherited from ``UsefulnessCriterionMatrix``

    #################### Trait metadata ####################
    ### ntrait                  inherited from ``UsefulnessCriterionMatrix``
    ### trait_axis              inherited from ``UsefulnessCriterionMatrix``

    ################ EPGC metadata property ################
    ### epgc                    inherited from ``UsefulnessCriterionMatrix``

    ########### Selection Percentile Properties ############
    ### upper_percentile        inherited from ``UsefulnessCriterionMatrix``
    ### selection_intensity     inherited from ``UsefulnessCriterionMatrix``

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    ### copy                    inherited from ``UsefulnessCriterionMatrix``
    ### deepcopy                inherited from ``UsefulnessCriterionMatrix``

    ######### Matrix element copy-on-manipulation ##########
    ### adjoin                  inherited from ``UsefulnessCriterionMatrix``
    ### delete                  inherited from ``UsefulnessCriterionMatrix``
    ### insert                  inherited from ``UsefulnessCriterionMatrix``
    ### select                  inherited from ``UsefulnessCriterionMatrix``
    ### adjoin_taxa             inherited from ``UsefulnessCriterionMatrix``
    ### delete_taxa             inherited from ``UsefulnessCriterionMatrix``
    ### insert_taxa             inherited from ``UsefulnessCriterionMatrix``
    ### select_taxa             inherited from ``UsefulnessCriterionMatrix``
    ### adjoin_trait            inherited from ``UsefulnessCriterionMatrix``
    ### delete_trait            inherited from ``UsefulnessCriterionMatrix``
    ### insert_trait            inherited from ``UsefulnessCriterionMatrix``
    ### select_trait            inherited from ``UsefulnessCriterionMatrix``

    ######### Matrix element in-place-manipulation #########
    ### append                  inherited from ``UsefulnessCriterionMatrix``
    ### remove                  inherited from ``UsefulnessCriterionMatrix``
    ### incorp                  inherited from ``UsefulnessCriterionMatrix``
    ### append_taxa             inherited from ``UsefulnessCriterionMatrix``
    ### remove_taxa             inherited from ``UsefulnessCriterionMatrix``
    ### incorp_taxa             inherited from ``UsefulnessCriterionMatrix``
    ### append_trait            inherited from ``UsefulnessCriterionMatrix``
    ### remove_trait            inherited from ``UsefulnessCriterionMatrix``
    ### incorp_trait            inherited from ``UsefulnessCriterionMatrix``

    ################### Sorting Methods ####################
    ### lexsort                 inherited from ``UsefulnessCriterionMatrix``
    ### reorder                 inherited from ``UsefulnessCriterionMatrix``
    ### sort                    inherited from ``UsefulnessCriterionMatrix``
    ### lexsort_taxa            inherited from ``UsefulnessCriterionMatrix``
    ### reorder_taxa            inherited from ``UsefulnessCriterionMatrix``
    ### sort_taxa               inherited from ``UsefulnessCriterionMatrix``
    ### lexsort_trait           inherited from ``UsefulnessCriterionMatrix``
    ### reorder_trait           inherited from ``UsefulnessCriterionMatrix``
    ### sort_trait              inherited from ``UsefulnessCriterionMatrix``

    ################### Grouping Methods ###################
    ### group                   inherited from ``UsefulnessCriterionMatrix``
    ### ungroup                 inherited from ``UsefulnessCriterionMatrix``
    ### is_grouped              inherited from ``UsefulnessCriterionMatrix``
    ### group_taxa              inherited from ``UsefulnessCriterionMatrix``
    ### ungroup_taxa            inherited from ``UsefulnessCriterionMatrix``
    ### is_grouped_taxa         inherited from ``UsefulnessCriterionMatrix``

    #################### Square Methods ####################
    ### is_square               inherited from ``UsefulnessCriterionMatrix``
    ### is_square_taxa          inherited from ``UsefulnessCriterionMatrix``

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat                  inherited from ``UsefulnessCriterionMatrix``
    ### concat_taxa             inherited from ``UsefulnessCriterionMatrix``
    ### concat_trait            inherited from ``UsefulnessCriterionMatrix``

    ################# Construction Methods #################
    ### from_gmod               inherited from ``UsefulnessCriterionMatrix``

    @classmethod
    @abstractmethod
    def from_algmod(
            cls, 
            algmod: AdditiveLinearGenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nmating: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            mem: int,
            **kwargs: dict
        ) -> 'AdditiveUsefulnessCriterionMatrix':
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
        mem : int
            Memory chunk size to use during matrix operations.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : AdditiveUsefulnessCriterionMatrix
            A matrix of additive usefulness criterion estimations.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_AdditiveUsefulnessCriterionMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``AdditiveUsefulnessCriterionMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, AdditiveUsefulnessCriterionMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                AdditiveUsefulnessCriterionMatrix.__name__,
                type(v).__name__
            )
        )
