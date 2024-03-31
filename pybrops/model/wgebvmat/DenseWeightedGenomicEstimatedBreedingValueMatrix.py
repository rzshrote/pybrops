"""
Module implementing dense wGEBV matrix representation.
"""

__all__ = [
    "DenseWeightedGenomicEstimatedBreedingValueMatrix",
    "check_is_DenseWeightedGenomicEstimatedBreedingValueMatrix",
]

from numbers import Real
from typing import Optional, Union
import numpy
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel, check_is_AdditiveLinearGenomicModel
from pybrops.model.wgebvmat.WeightedGenomicEstimatedBreedingValueMatrix import WeightedGenomicEstimatedBreedingValueMatrix
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix, check_is_GenotypeMatrix

class DenseWeightedGenomicEstimatedBreedingValueMatrix(
        DenseBreedingValueMatrix,
        WeightedGenomicEstimatedBreedingValueMatrix,
    ):
    """
    The ``DenseWeightedGenomicEstimatedBreedingValueMatrix`` class is used to represent 
    dense, multi-trait wGEBVs.

    Notes
    -----
    All elements within a ``DenseWeightedGenomicEstimatedBreedingValueMatrix``
    are mean-centered and scaled to unit variance for each trait.

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
        Constructor for the ``DenseWeightedGenomicEstimatedBreedingValueMatrix`` 
        class.
        
        Parameters
        ----------
        mat : numpy.ndarray
            An array of weighted genomic estimated breeding values of shape 
            ``(ntaxa,ntrait)``.

            Where:

            - ``ntaxa`` is the number of taxa.
            - ``ntrait`` is the number of traits.

            It is the responsibility of the user to ensure that the means and 
            standard deviations of this array along the ``taxa`` axis are ``0`` 
            and ``1``, respectively, if the breeding values are with respect to 
            the individuals in the breeding value matrix.

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
        # call DenseBreedingValueMatrix constructor
        super(DenseWeightedGenomicEstimatedBreedingValueMatrix, self).__init__(
            mat      = mat,
            location = location,
            scale    = scale,
            taxa     = taxa,
            taxa_grp = taxa_grp,
            trait    = trait,
            **kwargs
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ############################## Class Methods ###############################
    @classmethod
    def from_algmod(
            cls, 
            algmod: AdditiveLinearGenomicModel, 
            gmat: GenotypeMatrix, 
            **kwargs: dict
        ) -> 'DenseWeightedGenomicEstimatedBreedingValueMatrix':
        """
        Estimate Weighted Genomic Estimated Breeding Values (wGEBVs) from an 
        ``AdditiveLinearGenomicModel``.

        Parameters
        ----------
        algmod : AdditiveLinearGenomicModel
            An additive linear genomic model with which to estimate wGEBVs.
        
        gmat : GenotypeMatrix
            Genotypes for which to estimate wGEBVs.
        
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : DenseWeightedGenomicEstimatedBreedingValueMatrix
            A matrix object containing wGEBVs.
        """
        ### type checks
        check_is_AdditiveLinearGenomicModel(algmod, "algmod")
        check_is_GenotypeMatrix(gmat, "gmat")

        ### compute section
        # get genotype matrix as {0,1,2} format
        # (n,p)
        Z_a = gmat.mat_asformat("{0,1,2}")

        # get additive marker effects
        # (p,t)
        u_a = algmod.u_a

        # get favorable allele frequencies
        # (p,t)
        fafreq = algmod.fafreq(gmat)

        # get mask for where favorable allele frequency equals zero
        # (p,t)
        mask = fafreq == 0.0

        # calculate the numerator term: asin(1) - asin(sqrt(p))
        # scalar - (p,t) -> (p,t)
        numer = numpy.arcsin(1.0) - numpy.arcsin(numpy.sqrt(fafreq))

        # calculate p * (1-p)
        # (p,t)
        pq = fafreq * (1.0 - fafreq)

        # where p == 0, set to 1 to avoid division by zero
        # (p,t)
        pq[mask] = 1.0

        # calculate the denominator 1/sqrt(p * (1-p))
        # (p,t)
        denom = numpy.power(pq, -0.5)

        # calculate weights
        # (p,t) * (p,t) -> (p,t)
        weight = numer * denom

        # where p == 0, set to 1.0
        # (p,t)
        weight[mask] = 1.0

        # calculate wGEBVs
        # (p,t) * (p,t) -> (p,t)
        # (n,p) @ (p,t) -> (n,t)
        # (n,t)
        mat = Z_a @ (u_a * weight)

        # get metadata
        taxa = gmat.taxa
        taxa_grp = gmat.taxa_grp
        trait = algmod.trait

        # construct output
        out = cls.from_numpy(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
        )

        # attach metadata
        out.taxa_grp_name = gmat.taxa_grp_name
        out.taxa_grp_stix = gmat.taxa_grp_stix
        out.taxa_grp_spix = gmat.taxa_grp_spix
        out.taxa_grp_len  = gmat.taxa_grp_len

        return out

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_DenseWeightedGenomicEstimatedBreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``DenseWeightedGenomicEstimatedBreedingValueMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseWeightedGenomicEstimatedBreedingValueMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                DenseWeightedGenomicEstimatedBreedingValueMatrix.__name__,
                type(v).__name__
            )
        )
