"""
Module providing implementations of dense genotype matrices and associated error
checking routines.
"""

import copy
import numbers
from typing import Any, Optional, Sequence, Union
import cyvcf2
import h5py
import numpy
from numpy.typing import DTypeLike, ArrayLike

from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_io_h5py import check_group_in_hdf5
from pybrops.core.error.error_type_python import check_is_int
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_int8
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.DenseTaxaVariantMatrix import DenseTaxaVariantMatrix
from pybrops.core.util.h5py import save_dict_to_hdf5
from pybrops.popgen.gmap.DenseGeneticMappableMatrix import DenseGeneticMappableMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class DenseGenotypeMatrix(DenseTaxaVariantMatrix,DenseGeneticMappableMatrix,GenotypeMatrix):
    """
    A concrete class for unphased genoypte matrix objects.

    The purpose of this concrete class is to implement functionality for:
        1) Genotype matrix ploidy and phase metadata.
        2) Genotype matrix format conversion.
        3) Genotype matrix allele counting routines.
        4) Genotype matrix genotype counting routines.
        5) Loading genotype matrices from VCF and HDF5.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            vrnt_chrgrp: Optional[numpy.ndarray] = None,
            vrnt_phypos: Optional[numpy.ndarray] = None, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None,
            vrnt_xoprob: Optional[numpy.ndarray] = None, 
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None,
            vrnt_hapref: Optional[numpy.ndarray] = None, 
            vrnt_mask: Optional[numpy.ndarray] = None, 
            ploidy: int = 2, 
            **kwargs: dict
        ) -> None:
        """
        Parameters
        ----------
        mat : numpy.ndarray
            An int8 haplotype matrix. Must be {0,1,2} format.
        taxa : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.
        taxa_grp : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.
        vrnt_chrgrp : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant chromosome
            group labels. If ``None``, do not store any variant chromosome group
            label information.
        vrnt_phypos : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant chromosome
            physical positions. If ``None``, do not store any variant chromosome
            physical position information.
        vrnt_name : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant names.
            If ``None``, do not store any variant names.
        vrnt_genpos : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant chromosome
            genetic positions. If ``None``, do not store any variant chromosome
            genetic position information.
        vrnt_xoprob : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant crossover
            probabilities. If ``None``, do not store any variant crossover
            probabilities.
        vrnt_hapgrp : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant haplotype
            group labels. If ``None``, do not store any variant haplotype group
            label information.
        vrnt_hapalt : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant alternative
            alleles. If ``None``, do not store any variant alternative allele
            information.
        vrnt_hapref : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant reference
            alleles. If ``None``, do not store any variant reference allele
            information.
        vrnt_mask : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing a variant mask.
            If ``None``, do not store any variant mask information.
        ploidy : int
            The ploidy represented by the genotype matrix. This only represents
            ploidy of the reproductive habit. If the organism represented is an
            allopolyploid (e.g. hexaploid wheat), the ploidy is 2 since it
            reproduces in a diploid manner.
        """
        # set ploidy
        check_is_int(ploidy, "ploidy")
        self._ploidy = ploidy

        # TODO: Add a mat_format option to store as {0,1,2}, {-1,0,1}, etc.
        super(DenseGenotypeMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

    #################### Matrix copying ####################
    def __copy__(self):
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseGenotypeMatrix
            A shallow copy of the matrix.
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
            vrnt_chrgrp = copy.copy(self.vrnt_chrgrp),
            vrnt_phypos = copy.copy(self.vrnt_phypos),
            vrnt_name = copy.copy(self.vrnt_name),
            vrnt_genpos = copy.copy(self.vrnt_genpos),
            vrnt_xoprob = copy.copy(self.vrnt_xoprob),
            vrnt_hapgrp = copy.copy(self.vrnt_hapgrp),
            vrnt_hapalt = copy.copy(self.vrnt_hapalt),
            vrnt_hapref = copy.copy(self.vrnt_hapref),
            vrnt_mask = copy.copy(self.vrnt_mask),
            ploidy = self.ploidy
        )
        # copy taxa metadata
        out.taxa_grp_name = copy.copy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.copy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.copy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.copy(self.taxa_grp_len)

        # copy variant metadata
        out.vrnt_chrgrp_name = copy.copy(self.vrnt_chrgrp_name)
        out.vrnt_chrgrp_stix = copy.copy(self.vrnt_chrgrp_stix)
        out.vrnt_chrgrp_spix = copy.copy(self.vrnt_chrgrp_spix)
        out.vrnt_chrgrp_len = copy.copy(self.vrnt_chrgrp_len)

        return out

    def __deepcopy__(self, memo):
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict

        Returns
        -------
        out : Matrix
        """
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo),
            vrnt_chrgrp = copy.deepcopy(self.vrnt_chrgrp, memo),
            vrnt_phypos = copy.deepcopy(self.vrnt_phypos, memo),
            vrnt_name = copy.deepcopy(self.vrnt_name, memo),
            vrnt_genpos = copy.deepcopy(self.vrnt_genpos, memo),
            vrnt_xoprob = copy.deepcopy(self.vrnt_xoprob, memo),
            vrnt_hapgrp = copy.deepcopy(self.vrnt_hapgrp, memo),
            vrnt_hapalt = copy.deepcopy(self.vrnt_hapalt, memo),
            vrnt_hapref = copy.deepcopy(self.vrnt_hapref, memo),
            vrnt_mask = copy.deepcopy(self.vrnt_mask, memo),
            ploidy = self.ploidy
        )
        # copy taxa metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name, memo)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix, memo)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix, memo)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len, memo)

        # copy variant metadata
        out.vrnt_chrgrp_name = copy.deepcopy(self.vrnt_chrgrp_name, memo)
        out.vrnt_chrgrp_stix = copy.deepcopy(self.vrnt_chrgrp_stix, memo)
        out.vrnt_chrgrp_spix = copy.deepcopy(self.vrnt_chrgrp_spix, memo)
        out.vrnt_chrgrp_len = copy.deepcopy(self.vrnt_chrgrp_len, memo)

        return out

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Genotype Data Properites ##############
    @DenseTaxaVariantMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        check_is_ndarray(value, "mat")
        check_ndarray_dtype_is_int8(value, "mat")
        check_ndarray_ndim(value, "mat", 2)
        # check_ndarray_in_interval(value, "mat", 0, self._ploidy + 1)
        self._mat = value

    ############## General matrix properties ###############
    @property
    def ploidy(self) -> int:
        """Genome ploidy number represented by matrix."""
        return self._ploidy
    @ploidy.setter
    def ploidy(self, value: int) -> None:
        """Set matrix ploidy number"""
        error_readonly("ploidy")
    @ploidy.deleter
    def ploidy(self) -> None:
        """Delete matrix ploidy number"""
        error_readonly("ploidy")

    @property
    def nphase(self) -> int:
        """Number of chromosome phases represented by the matrix."""
        return 0
    @nphase.setter
    def nphase(self, value: int) -> None:
        """Set number of phases"""
        error_readonly("nphase")
    @nphase.deleter
    def nphase(self) -> None:
        """Delete number of phases"""
        error_readonly("nphase")

    @property
    def mat_format(self) -> str:
        """Matrix representation format."""
        return "{0,1,2}"
    @mat_format.setter
    def mat_format(self, value: str) -> None:
        """Set matrix representation format"""
        error_readonly("mat_format")
    @mat_format.deleter
    def mat_format(self) -> None:
        """Delete matrix representation format"""
        error_readonly("mat_format")

    ############### Taxa Metadata Properites ###############
    @DenseTaxaVariantMatrix.taxa_axis.getter
    def taxa_axis(self) -> int:
        """Get taxa axis number"""
        return 0

    ############# Variant Metadata Properites ##############
    @DenseTaxaVariantMatrix.vrnt_axis.getter
    def vrnt_axis(self) -> int:
        """Get variant axis"""
        return 1

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin_taxa(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseGenotypeMatrix':
        """
        Add additional elements to the end of the Matrix along the taxa axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        out = super(DenseGenotypeMatrix, self).adjoin_taxa(
            values = values,
            taxa = taxa,
            taxa_grp = taxa_grp,
            ploidy = self.ploidy,
            **kwargs
        )

        return out

    def adjoin_vrnt(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            vrnt_chrgrp: Optional[numpy.ndarray] = None, 
            vrnt_phypos: Optional[numpy.ndarray] = None, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None, 
            vrnt_xoprob: Optional[numpy.ndarray] = None, 
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None, 
            vrnt_hapref: Optional[numpy.ndarray] = None, 
            vrnt_mask: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseGenotypeMatrix':
        """
        Add additional elements to the end of the Matrix along the variant axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to adjoin to the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to adjoin to the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to adjoin to the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to adjoin to the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to adjoin to the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to adjoin to the Matrix.
        vrnt_hapalt : numpy.ndarray
            Variant haplotype sequence.
        vrnt_hapref : numpy.ndarray
            Variant haplotype reference sequence.
        vrnt_mask : numpy.ndarray
            Variant mask to adjoin to the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : Matrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        out = super(DenseGenotypeMatrix, self).adjoin_vrnt(
            values = values,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
            ploidy = self.ploidy,
            **kwargs
        )

        return out

    def delete_taxa(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'DenseGenotypeMatrix':
        """
        Delete sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseGenotypeMatrix
            A DenseGenotypeMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseGenotypeMatrix is allocated and filled.
        """
        out = super(DenseGenotypeMatrix, self).delete_taxa(
            obj = obj,
            ploidy = self.ploidy,
            **kwargs
        )

        return out

    def delete_vrnt(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'DenseGenotypeMatrix':
        """
        Delete sub-arrays along the variant axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseGenotypeMatrix
            A DenseGenotypeMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseGenotypeMatrix is allocated and filled.
        """
        out = super(DenseGenotypeMatrix, self).delete_vrnt(
            obj = obj,
            ploidy = self.ploidy,
            **kwargs
        )

        return out

    def insert_taxa(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseGenotypeMatrix':
        """
        Insert values along the taxa axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        taxa : numpy.ndarray
            Taxa names to insert into the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseGenotypeMatrix
            A DenseGenotypeMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseGenotypeMatrix is allocated and filled.
        """
        # create output
        out = super(DenseGenotypeMatrix, self).insert_taxa(
            obj = obj,
            values = values,
            taxa = taxa,
            taxa_grp = taxa_grp,
            ploidy = self.ploidy,
            **kwargs
        )

        return out

    def insert_vrnt(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            vrnt_chrgrp: Optional[numpy.ndarray] = None, 
            vrnt_phypos: Optional[numpy.ndarray] = None, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None, 
            vrnt_xoprob: Optional[numpy.ndarray] = None, 
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None, 
            vrnt_hapref: Optional[numpy.ndarray] = None, 
            vrnt_mask: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseGenotypeMatrix':
        """
        Insert values along the variant axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : array_like
            Values to insert into the matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to insert into the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to insert into the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to insert into the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to insert into the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to insert into the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to insert into the Matrix.
        vrnt_mask : numpy.ndarray
            Variant mask to insert into the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseGenotypeMatrix
            A DenseGenotypeMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseGenotypeMatrix is allocated and filled.
        """
        # create output
        out = super(DenseGenotypeMatrix, self).insert_vrnt(
            obj = obj,
            values = values,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
            ploidy = self.ploidy,
            **kwargs
        )

        return out

    def select_taxa(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'DenseGenotypeMatrix':
        """
        Select certain values from the Matrix along the taxa axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseGenotypeMatrix
            The output DenseGenotypeMatrix with values selected. Note that select does not
            occur in-place: a new DenseGenotypeMatrix is allocated and filled.
        """
        out = super(DenseGenotypeMatrix, self).select_taxa(
            indices = indices,
            ploidy = self.ploidy,
            **kwargs
        )

        return out

    def select_vrnt(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'DenseGenotypeMatrix':
        """
        Select certain values from the Matrix along the variant axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseGenotypeMatrix
            The output DenseGenotypeMatrix with values selected. Note that select does not
            occur in-place: a new DenseGenotypeMatrix is allocated and filled.
        """
        out = super(DenseGenotypeMatrix, self).select_vrnt(
            indices = indices,
            ploidy = self.ploidy,
            **kwargs
        )

        return out

    @classmethod
    def concat_taxa(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DenseGenotypeMatrix':
        """
        Concatenate list of Matrix together along the taxa axis.

        Parameters
        ----------
        mats : Sequence of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : DenseGenotypeMatrix
            The concatenated DenseGenotypeMatrix. Note that concat does not occur in-place:
            a new DenseGenotypeMatrix is allocated and filled.
        """
        out = super(DenseGenotypeMatrix, cls).concat_taxa(
            mats = mats,
            ploidy = mats[0].ploidy,
            **kwargs
        )

        # copy metadata from source
        out.vrnt_chrgrp_name = mats[0].vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = mats[0].vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = mats[0].vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = mats[0].vrnt_chrgrp_len

        return out

    @classmethod
    def concat_vrnt(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DenseGenotypeMatrix':
        """
        Concatenate list of Matrix together along the variant axis.

        Parameters
        ----------
        mats : Sequence of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        out = super(DenseGenotypeMatrix, cls).concat_vrnt(
            mats = mats,
            ploidy = mats[0].ploidy,
            **kwargs
        )

        return out

    ################## Matrix conversion ###################
    def mat_asformat(
            self, 
            format: str
        ) -> numpy.ndarray:
        """
        Get mat in a specific format type.

        Parameters
        ---------
        format : str
            Desired output format. Options are "{0,1,2}", "{-1,0,1}", "{-1,m,1}".

        Returns
        -------
        out : numpy.ndarray
            Matrix in the desired output format.
        """
        if format == "{0,1,2}":
            return self.mat.copy()
        elif format == "{-1,0,1}":
            out = self.mat - 1
            return out
        elif format == "{-1,m,1}":
            # OPTIMIZE: there's probably a matrix multiplication way to do this instead of if-else
            out = self.mat - 1.0    # (n,p) float64
            for i in range(out.shape[1]):
                view = out[:,i]
                mean = view.mean()
                mask = (view == 0)
                out[mask,i] = mean
            return out
        else:
            raise ValueError('Format not recognized. Options are "{0,1,2}", "{-1,0,1}", "{-1,m,1}".')

    ############## Matrix summary statistics ###############
    def tacount(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Allele count of the non-zero allele within each taxon.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the accumulator and returned array.
            If ``None``, use the native accumulator type (int or float).

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(n,p)`` containing allele counts of the
            allele coded as ``1`` for all ``n`` individuals, for all ``p`` loci.
        """
        # get accumulator type
        if dtype is None:
            if numpy.issubdtype(self._mat.dtype, numpy.integer):
                dtype = int
            elif numpy.issubdtype(self._mat.dtype, numpy.floating):
                dtype = float
            else:
                raise ValueError("No default accumulator type for GenotypeMatrix dtype {0}".format(self._mat.dtype))
        # calculate the taxa allele sums
        out = self._mat.astype(dtype)
        return out

    def tafreq(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Allele frequency of the non-zero allele within each taxon.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(n,p)`` containing allele frequencies of
            the allele coded as ``1`` for all ``n`` individuals, for all ``p``
            loci.
        """
        rnphase = 1.0 / self.ploidy         # take the reciprocal of the ploidy number
        out = rnphase * self._mat           # take sum across the phase axis (0) and divide by ploidy
        if dtype is not None:               # if dtype is specified
            dtype = numpy.dtype(dtype)      # ensure conversion to dtype class
            if out.dtype != dtype:          # if output dtype and desired are different
                out = dtype.type(out)       # convert to correct dtype
        return out

    def acount(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Allele count of the non-zero allele across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the accumulator and returned array. If ``None``, use the native integer type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing allele counts of the
            allele coded as ``1`` for all ``p`` loci.
        """
        # process dtype
        if dtype is None:
            dtype = int
        dtype = numpy.dtype(dtype)

        # take sum across the taxa axis
        out = self._mat.sum(self.taxa_axis, dtype = dtype)
        
        return out

    def afreq(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Allele frequency of the non-zero allele across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing allele frequencies of
            the allele coded as ``1`` for all ``p`` loci.
        """
        denom = (self.ploidy * self.ntaxa)              # get ploidy * ntaxa
        rnphase = 1.0 / denom                           # take 1 / (ploidy * ntaxa)
        out = rnphase * self._mat.sum(self.taxa_axis)   # take sum across the taxa axis (0) and divide by nphase
        if dtype is not None:                           # if dtype is specified
            dtype = numpy.dtype(dtype)                  # ensure conversion to dtype class
            if out.dtype != dtype:                      # if output dtype and desired are different
                out = dtype.type(out)                   # convert to correct dtype
        return out

    def apoly(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Allele polymorphism presence or absense across all loci.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing indicator variables for
            whether the locus is polymorphic.
        """
        # find any loci that are monomorphic for 0 allele
        min_mask = numpy.all(self.mat == 0, axis = self.taxa_axis)

        # find any loci that are monomorphic for 1 allele
        max_mask = numpy.all(self.mat == self.ploidy, axis = self.taxa_axis)

        # logical or results and take logical not
        out = numpy.logical_not(min_mask | max_mask)

        # convert to specific dtype if needed
        if dtype is not None:                           # if dtype is specified
            dtype = numpy.dtype(dtype)                  # ensure conversion to dtype class
            if out.dtype != dtype:                      # if output dtype and desired are different
                out = dtype.type(out)                   # convert to correct dtype

        return out

    def maf(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Minor allele frequency across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing allele frequencies for
            the minor allele.
        """
        out = self.afreq(dtype)     # get allele frequencies
        mask = out > 0.5            # create mask of allele frequencies > 0.5
        out[mask] = 1.0 - out[mask] # take 1 - allele frequency
        return out

    def meh(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numbers.Number:
        """
        Mean expected heterozygosity across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numbers.Number
            A number representing the mean expected heterozygous.
            If ``dtype`` is ``None``, then a native 64-bit floating point is
            returned. Otherwise, of type specified by ``dtype``.
        """
        # OPTIMIZE: take dot product with allele counts, then divide?
        p = self.afreq()                    # get haplotype frequency (p)
        out = numpy.dot(p, 1.0 - p)         # take p*(1-p) across all loci and sum the products
        rnphase = self.ploidy / self.nvrnt  # ploidy / nvrnt
        out *= rnphase                      # multiply summation by (nphase / nvrnt)
        if dtype is not None:               # if dtype is specified
            dtype = numpy.dtype(dtype)      # ensure conversion to dtype class
            if out.dtype != dtype:          # if output dtype and desired are different
                out = dtype.type(out)       # convert to correct dtype
        return out

    def gtcount(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Gather genotype counts for homozygous major, heterozygous, homozygous
        minor for all individuals.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray array of shape ``(g,p)`` containing allele counts
            across all ``p`` loci for each of ``g`` genotype combinations.

            Where:

            - ``out[0]`` is the count of ``0`` genotype across all loci
            - ``out[1]`` is the count of ``1`` genotype across all loci
            - ``out[2]`` is the count of ``2`` genotype across all loci
            - ``...``
            - ``out[g-1]`` is the count of ``g-1`` genotype across all loci
        """
        ngt = self.nphase + 1   # get number of genotype combos
        mat = self._mat         # get matrix

        out = numpy.empty(      # allocate output array
            (ngt, self.nvrnt),  # shape (g, p)
            dtype='int64'       # int64 dtype
        )

        # get correct axis to sum across
        axis = self.taxa_axis

        for i in range(ngt):                # for each genotype combo
            out[i] = (mat == i).sum(axis)   # record counts for genotype 'i'

        if dtype is not None:               # if dtype is specified
            dtype = numpy.dtype(dtype)      # ensure conversion to dtype class
            if out.dtype != dtype:          # if output dtype and desired are different
                out = dtype.type(out)       # convert to correct dtype

        return out

    def gtfreq(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Gather genotype frequencies for homozygous major, heterozygous,
        homozygous minor across all individuals.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray array of shape ``(g,p)`` containing haplotype counts
            across all ``p`` loci for all ``g`` genotype combinations.

            Where:

            - ``out[0]`` is the frequency of ``0`` genotype across all loci
            - ``out[1]`` is the frequency of ``1`` genotype across all loci
            - ``out[2]`` is the frequency of ``2`` genotype across all loci
            - ``...``
            - ``out[g-1]`` is the frequency of ``g-1`` genotype across all loci
        """
        recip = 1.0 / self.ntaxa        # get reciprocal of number of taxa
        out = recip * self.gtcount()    # calculate genotype frequencies
        if dtype is not None:           # if dtype is specified
            dtype = numpy.dtype(dtype)  # ensure conversion to dtype class
            if out.dtype != dtype:      # if output dtype and desired are different
                out = dtype.type(out)   # convert to correct dtype
        return out

    ################### Matrix File I/O ####################
    def to_hdf5(
            self, 
            filename: str, 
            groupname: Optional[str] = None
        ) -> None:
        """
        Write GenotypeMatrix to an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name to which to write.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is written to the base HDF5 group.
        """
        h5file = h5py.File(filename, "a")                       # open HDF5 in write mode
        ######################################################### process groupname argument
        if isinstance(groupname, str):                          # if we have a string
            if groupname[-1] != '/':                            # if last character in string is not '/'
                groupname += '/'                                # add '/' to end of string
        elif groupname is None:                                 # else if groupname is None
            groupname = ""                                      # empty string
        else:                                                   # else raise error
            raise TypeError("'groupname' must be of type str or None")
        ######################################################### populate HDF5 file
        data_dict = {                                           # data dictionary
            "mat": self.mat,
            "taxa" : self.taxa,
            "taxa_grp" : self.taxa_grp,
            "vrnt_chrgrp" : self.vrnt_chrgrp,
            "vrnt_phypos" : self.vrnt_phypos,
            "vrnt_name" : self.vrnt_name,
            "vrnt_genpos" : self.vrnt_genpos,
            "vrnt_xoprob" : self.vrnt_xoprob,
            "vrnt_hapgrp" : self.vrnt_hapgrp,
            "vrnt_hapalt" : self.vrnt_hapalt,
            "vrnt_hapref" : self.vrnt_hapref,
            "vrnt_mask" : self.vrnt_mask,
            "ploidy" : self.ploidy
        }
        save_dict_to_hdf5(h5file, groupname, data_dict)         # save data
        ######################################################### write conclusion
        h5file.close()                                          # close the file

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(
            cls, 
            filename: str, 
            groupname: Optional[str] = None
        ) -> 'DenseGenotypeMatrix':
        """
        Read GenotypeMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is read from base HDF5 group.

        Returns
        -------
        gmat : GenotypeMatrix
            A genotype matrix read from file.
        """
        check_file_exists(filename)                             # check file exists
        h5file = h5py.File(filename, "r")                       # open HDF5 in read only
        ######################################################### process groupname argument
        if isinstance(groupname, str):                          # if we have a string
            check_group_in_hdf5(groupname, h5file, filename)    # check that group exists
            if groupname[-1] != '/':                            # if last character in string is not '/'
                groupname += '/'                                # add '/' to end of string
        elif groupname is None:                                 # else if groupname is None
            groupname = ""                                      # empty string
        else:                                                   # else raise error
            raise TypeError("'groupname' must be of type str or None")
        ######################################################### check that we have all required fields
        required_fields = ["mat"]                               # all required arguments
        for field in required_fields:                           # for each required field
            fieldname = groupname + field                       # concatenate base groupname and field
            check_group_in_hdf5(fieldname, h5file, filename)    # check that group exists
        ######################################################### read data
        data_dict = {                                           # output dictionary
            "mat": None,
            "taxa" : None,
            "taxa_grp" : None,
            "vrnt_chrgrp" : None,
            "vrnt_phypos" : None,
            "vrnt_name" : None,
            "vrnt_genpos" : None,
            "vrnt_xoprob" : None,
            "vrnt_hapgrp" : None,
            "vrnt_hapalt" : None,
            "vrnt_hapref" : None,
            "vrnt_mask" : None,
            "ploidy" : None
        }
        for field in data_dict.keys():                          # for each field
            fieldname = groupname + field                       # concatenate base groupname and field
            if fieldname in h5file:                             # if the field exists in the HDF5 file
                data_dict[field] = h5file[fieldname][()]        # read array
        # fieldname = groupname + "ploidy"                        # special case for reading ploidy ()
        # if fieldname in h5file:
        #     data_dict["ploidy"] = h5file[fieldname][]
        ######################################################### read conclusion
        h5file.close()                                          # close file
        if data_dict["taxa"] is not None:
            data_dict["taxa"] = numpy.object_(                  # convert taxa strings from byte to utf-8
                [s.decode("utf-8") for s in data_dict["taxa"]]
            )
        if data_dict["vrnt_name"] is not None:
            data_dict["vrnt_name"] = numpy.object_(             # convert vrnt_name string from byte to utf-8
                [s.decode("utf-8") for s in data_dict["vrnt_name"]]
            )
        if data_dict["vrnt_hapgrp"] is not None:
            data_dict["vrnt_hapgrp"] = numpy.object_(           # convert vrnt_hapgrp string from byte to utf-8
                [s.decode("utf-8") for s in data_dict["vrnt_hapgrp"]]
            )
        if data_dict["vrnt_hapalt"] is not None:
            data_dict["vrnt_hapalt"] = numpy.object_(           # convert vrnt_hapalt string from byte to utf-8
                [s.decode("utf-8") for s in data_dict["vrnt_hapalt"]]
            )
        if data_dict["ploidy"] is not None:
            data_dict["ploidy"] = int(data_dict["ploidy"])      # convert to int, from numpy.int64
        ######################################################### create object
        gmat = cls(**data_dict)                                 # create object from read data
        return gmat

    @classmethod
    def from_vcf(
            cls, 
            fname: str
        ) -> 'DenseGenotypeMatrix':
        """
        Create a DenseGenotypeMatrix from a VCF file.

        Parameters
        ----------
        fname : str
            Path to VCF file.

        Returns
        -------
        out : DenseGenotypeMatrix
            An unphased genotype matrix with associated metadata from VCF file.
        """
        # make VCF iterator
        vcf = cyvcf2.VCF(fname)

        # extract taxa names from vcf header
        taxa = numpy.object_(vcf.samples)

        # make empty lists to store extracted values
        mat = []
        vrnt_chrgrp = []
        vrnt_phypos = []
        vrnt_name = []

        # iterate through VCF file and accumulate variants
        for variant in vcf:
            # append chromosome integer
            vrnt_chrgrp.append(int(variant.CHROM))

            # append variant position coordinates
            vrnt_phypos.append(variant.POS)

            # append marker name
            vrnt_name.append(str(variant.ID))

            # extract allele states + whether they are phased or not
            phases = numpy.int8(variant.genotypes)

            # check that they are all phased
            #pybrops.util.check_matrix_all_value(phases[:,2], "is_phased", True)

            # TODO: maybe modify shapes here to avoid transpose and copy below?
            # append genotype states
            mat.append(phases[:,0:2].copy())

        # convert to compact int8 type and transpose genotype matrix
        # (p,n,m) --transpose-> (m,n,p)
        mat = numpy.int8(mat).transpose(2,1,0)

        # extract ploidy of sample
        ploidy = mat.shape[0]

        # sum across phases
        # (m,n,p).sum(0) -> (n,p)
        mat = mat.sum(0, dtype = 'int8')

        # convert to numpy.ndarray
        vrnt_chrgrp = numpy.int64(vrnt_chrgrp)  # convert to int64 array
        vrnt_phypos = numpy.int64(vrnt_phypos)  # convert to int64 array
        vrnt_name = numpy.object_(vrnt_name)    # convert to object array

        pvm = cls(
            mat = mat,
            taxa = taxa,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            ploidy = ploidy
        )

        return pvm



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseGenotypeMatrix(v: object) -> bool:
    """
    Determine whether an object is a ``DenseGenotypeMatrix``.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``v`` is a ``DenseGenotypeMatrix`` object instance.
    """
    return isinstance(v, DenseGenotypeMatrix)

def check_is_DenseGenotypeMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``DenseGenotypeMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseGenotypeMatrix):
        raise TypeError("'{0}' must be a DenseGenotypeMatrix.".format(varname))
