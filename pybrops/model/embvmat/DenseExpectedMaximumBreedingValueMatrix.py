"""
Module implementing dense EMBV matrix representation.
"""

__all__ = [
    "DenseExpectedMaximumBreedingValueMatrix",
    "check_is_DenseExpectedMaximumBreedingValueMatrix",
]

from numbers import Integral, Real
from typing import Optional, Union
import numpy
from pybrops.core.util.mate import dense_dh
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_shape_eq
from pybrops.core.random.prng import global_prng
from pybrops.model.embvmat.ExpectedMaximumBreedingValueMatrix import ExpectedMaximumBreedingValueMatrix
from pybrops.model.gmod.GenomicModel import GenomicModel, check_is_GenomicModel
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.popgen.gmat.DensePhasedGenotypeMatrix import DensePhasedGenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_PhasedGenotypeMatrix_has_vrnt_xoprob, check_is_PhasedGenotypeMatrix

class DenseExpectedMaximumBreedingValueMatrix(
        DenseBreedingValueMatrix,
        ExpectedMaximumBreedingValueMatrix
    ):
    """
    The ``DenseExpectedMaximumBreedingValueMatrix`` class is used to represent 
    dense, multi-trait EMBVs.

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
        Constructor for DenseExpectedMaximumBreedingValueMatrix.
        
        Parameters
        ----------
        mat : numpy.ndarray
            An array of expected maximum breeding values of shape ``(ntaxa,ntrait)``.

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
        # call DenseBreedingValueMatrix constructor
        super(DenseExpectedMaximumBreedingValueMatrix, self).__init__(
            mat = mat,
            location = location,
            scale = scale,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ############################## Class Methods ###############################
    @classmethod
    def from_gmod(
            cls,
            gmod: GenomicModel,
            pgmat: PhasedGenotypeMatrix,
            nprogeny: Union[Integral,numpy.ndarray],
            nrep: Union[Integral,numpy.ndarray],
            **kwargs: dict
        ) -> 'ExpectedMaximumBreedingValueMatrix':
        """
        Estimate Expected Maximum Breeding Values (EMBVs) from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            A ``GenomicModel`` with which to estimate EMBVs.

        pgmat : PhasedGenotypeMatrix
            Genomes for which to estimate EMBVs.

        nprogeny : Integral, numpy.ndarray
            Number of DH progenies to create per individual.
            If ``Integral``, each taxon produces ``nprogeny`` DH progenies.
            If ``numpy.ndarray`` of shape ``(ntaxa,)``, each taxon produces the
            number of progenies corresponding to elements in the array.

        nrep : Integral, numpy.ndarray
            Number of sample replicates for which to sample ``nprogeny`` DH progenies.
            If ``Integral``, each taxon produces ``nrep`` sample replicates.
            If ``numpy.ndarray`` of shape ``(ntaxa,)``, each taxon produces the
            number of sample replicates corresponding to elements in the array.

        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : ExpectedMaximumBreedingValueMatrix
            An ``ExpectedMaximumBreedingValueMatrix`` object containing EMBVs.
        """
        ################### check inputs ###################
        # check gmod
        check_is_GenomicModel(gmod, "gmod")

        # check pgmat
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_PhasedGenotypeMatrix_has_vrnt_xoprob(pgmat, "pgmat")

        # check nprogeny
        if isinstance(nprogeny, Integral):
            nprogeny = numpy.repeat(nprogeny, pgmat.ntaxa)
        check_is_ndarray(nprogeny, "nprogeny")
        check_ndarray_shape_eq(nprogeny, "nprogeny", (pgmat.ntaxa,))

        # check nrep
        if isinstance(nrep, Integral):
            nrep = numpy.repeat(nrep, pgmat.ntaxa)
        check_is_ndarray(nrep, "nrep")
        check_ndarray_shape_eq(nrep, "nrep", (pgmat.ntaxa,))

        ################# calculate EMBVs ##################
        # get pointers to important arrays
        geno        = pgmat.mat
        vrnt_chrgrp = pgmat.vrnt_chrgrp
        vrnt_phypos = pgmat.vrnt_phypos
        vrnt_name   = pgmat.vrnt_name
        vrnt_genpos = pgmat.vrnt_genpos
        vrnt_xoprob = pgmat.vrnt_xoprob
        vrnt_hapgrp = pgmat.vrnt_hapgrp
        vrnt_hapalt = pgmat.vrnt_hapalt
        vrnt_hapref = pgmat.vrnt_hapref
        vrnt_mask   = pgmat.vrnt_mask

        # allocate an array to store EMBVs
        embv = numpy.empty((pgmat.ntaxa,gmod.ntrait), dtype = float)

        # for each taxon
        for i in range(pgmat.ntaxa):

            # allocate an array to store maximum breeding values
            # (nrep,ntrait)
            mbv = numpy.empty((nrep[i],gmod.ntrait))

            # for each replicate
            for j in range(nrep[i]):
            
                # create dh progenies
                mat = dense_dh(geno, numpy.repeat(i, nprogeny[i]), vrnt_xoprob, global_prng)

                # TODO: remove dependency on implementation
                # create genotype matrix
                progeny = DensePhasedGenotypeMatrix(
                    mat         = mat,
                    vrnt_chrgrp = vrnt_chrgrp,
                    vrnt_phypos = vrnt_phypos,
                    vrnt_name   = vrnt_name,
                    vrnt_genpos = vrnt_genpos,
                    vrnt_xoprob = vrnt_xoprob,
                    vrnt_hapgrp = vrnt_hapgrp,
                    vrnt_hapalt = vrnt_hapalt,
                    vrnt_hapref = vrnt_hapref,
                    vrnt_mask   = vrnt_mask,
                )

                # predict progeny breeding values
                # (nprogeny,t)
                bvmat = gmod.gebv(progeny)

                # store maximum trait values
                mbv[j,:] = bvmat.tmax(unscale = True)
            
            # take the mean maximum breeding value and stor into embv
            # (i,ntrait) <- (ntrait,)
            embv[i,:] = mbv.mean(axis = 0)

        # construct a breeding value matrix
        out = cls.from_numpy(
            mat = embv,
            taxa = pgmat.taxa,
            taxa_grp = pgmat.taxa_grp,
            trait = gmod.trait,
        )

        return out

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_DenseExpectedMaximumBreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``DenseExpectedMaximumBreedingValueMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseExpectedMaximumBreedingValueMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                DenseExpectedMaximumBreedingValueMatrix.__name__,
                type(v).__name__
            )
        )
