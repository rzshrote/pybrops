"""
Module providing implementations of dense genotype matrices and associated error
checking routines.
"""

__all__ = [
    "DenseGenotypeMatrix",
    "check_is_DenseGenotypeMatrix",
]

import copy
from numbers import Real
from pathlib import Path
from typing import Optional
from typing import Sequence
from typing import Union
import cyvcf2
import h5py
import numpy
from numpy.typing import DTypeLike
from numpy.typing import ArrayLike

from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_type_python import check_is_int
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_int8
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group
from pybrops.core.error.error_value_h5py import check_h5py_File_is_readable
from pybrops.core.error.error_value_h5py import check_h5py_File_is_writable
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.DenseTaxaVariantMatrix import DenseTaxaVariantMatrix
from pybrops.core.util.h5py import h5py_File_read_int
from pybrops.core.util.h5py import h5py_File_read_ndarray
from pybrops.core.util.h5py import h5py_File_read_ndarray_int8
from pybrops.core.util.h5py import h5py_File_read_ndarray_utf8
from pybrops.core.util.h5py import h5py_File_write_dict
from pybrops.popgen.gmap.DenseGeneticMappableMatrix import DenseGeneticMappableMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class DenseGenotypeMatrix(
        DenseTaxaVariantMatrix,
        DenseGeneticMappableMatrix,
        GenotypeMatrix,
    ):
    """
    A concrete class for unphased genoypte matrix objects.

    The purpose of this concrete class is to implement functionality for:
        1) Genotype matrix ploidy and phase metadata.
        2) Genotype matrix format conversion.
        3) Genotype matrix allele counting routines.
        4) Genotype matrix genotype counting routines.
        5) Loading genotype matrices from VCF and HDF5.
    """

    ########################## Special Object Methods ##########################
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
            A ``numpy.ndarray`` of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.
        
            Where:

            - ``n`` is the number of taxa.

        taxa_grp : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.
        
            Where:

            - ``n`` is the number of taxa.

        vrnt_chrgrp : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(p,)`` containing variant chromosome
            group labels. If ``None``, do not store any variant chromosome group
            label information.

            Where:

            - ``p`` is the number of variants.

        vrnt_phypos : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(p,)`` containing variant chromosome
            physical positions. If ``None``, do not store any variant chromosome
            physical position information.
        
            Where:

            - ``p`` is the number of variants.

        vrnt_name : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(p,)`` containing variant names.
            If ``None``, do not store any variant names.
        
            Where:

            - ``p`` is the number of variants.

        vrnt_genpos : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(p,)`` containing variant chromosome
            genetic positions. If ``None``, do not store any variant chromosome
            genetic position information.
        
            Where:

            - ``p`` is the number of variants.

        vrnt_xoprob : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(p,)`` containing variant crossover
            probabilities. If ``None``, do not store any variant crossover
            probabilities.
        
            Where:

            - ``p`` is the number of variants.

        vrnt_hapgrp : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(p,)`` containing variant haplotype
            group labels. If ``None``, do not store any variant haplotype group
            label information.
        
            Where:

            - ``p`` is the number of variants.

        vrnt_hapalt : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(p,)`` containing variant alternative
            alleles. If ``None``, do not store any variant alternative allele
            information.
        
            Where:

            - ``p`` is the number of variants.

        vrnt_hapref : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(p,)`` containing variant reference
            alleles. If ``None``, do not store any variant reference allele
            information.
        
            Where:

            - ``p`` is the number of variants.

        vrnt_mask : numpy.ndarray, None
            A ``numpy.ndarray`` of shape ``(p,)`` containing a variant mask.
            If ``None``, do not store any variant mask information.
        
            Where:

            - ``p`` is the number of variants.

        ploidy : int
            The ploidy represented by the genotype matrix. This only represents
            ploidy of the reproductive habit. If the organism represented is an
            allopolyploid (e.g. hexaploid wheat), the ploidy is 2 since it
            reproduces in a diploid manner.
        
        kwargs : dict
            Additional keyword arguments.
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
    def __copy__(
            self
        ) -> 'DenseGenotypeMatrix':
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

    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'DenseGenotypeMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

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
        return "<{0} of shape (nphase = {1}, ntaxa = {2}, nvrnt = {3}) with ploidy = {4} at {5}>".format(
            type(self).__name__,
            self.nphase,
            self.ntaxa,
            self.nvrnt,
            self.ploidy,
            hex(id(self)),
        )

    ############################ Object Properties #############################

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

    @property
    def nphase(self) -> int:
        """Number of chromosome phases represented by the matrix."""
        return 0
    @nphase.setter
    def nphase(self, value: int) -> None:
        """Set number of phases"""
        error_readonly("nphase")

    @property
    def mat_format(self) -> str:
        """Matrix representation format."""
        return "{0,1,2}"
    @mat_format.setter
    def mat_format(self, value: str) -> None:
        """Set matrix representation format"""
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

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseGenotypeMatrix':
        """
        Make a shallow copy of the DenseGenotypeMatrix.

        Returns
        -------
        out : DenseGenotypeMatrix
            A shallow copy of the original DenseGenotypeMatrix.
        """
        return self.__copy__()

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseGenotypeMatrix':
        """
        Make a deep copy of the DenseGenotypeMatrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseGenotypeMatrix
            A deep copy of the original DenseGenotypeMatrix.
        """
        return self.__deepcopy__(memo)

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
        vrnt_hapalt : numpy.ndarray
            Variant alternative haplotype labels to insert into the Matrix.
        vrnt_hapref : numpy.ndarray
            Variant reference haplotype labels to insert into the Matrix.
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
            The data type of the accumulator and returned array.
            If ``None``, use the native accumulator type (int or float).

        Returns
        -------
        out : numpy.ndarray
            A ``numpy.ndarray`` of shape ``(n,p)`` containing allele counts of the
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
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A ``numpy.ndarray`` of shape ``(n,p)`` containing allele frequencies of
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
            The data type of the accumulator and returned array. If ``None``, use the native integer type.

        Returns
        -------
        out : numpy.ndarray
            A ``numpy.ndarray`` of shape ``(p,)`` containing allele counts of the
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
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A ``numpy.ndarray`` of shape ``(p,)`` containing allele frequencies of
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

    def afixed(
            self,
            dtype: Optional[DTypeLike] = None,
        ) -> numpy.ndarray:
        """
        Determine allele fixation for loci across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The data type of the returned array. 
            If ``None``, use the native type.
        
        Returns
        -------
        out : numpy.ndarray
            A ``numpy.ndarray`` of shape ``(p,)`` containing indicator variables 
            for whether a locus is fixed at a particular locus.
        """
        # get the allele frequency
        afreq = self.afreq()

        # determine whether the allele frequency is equal to 0.0 or 1.0
        out = (afreq == 0.0) | (afreq == 1.0)

        # convert to specific dtype if needed
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
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A ``numpy.ndarray`` of shape ``(p,)`` containing indicator variables for
            whether the locus is polymorphic.
        """
        # get the allele frequency
        afreq = self.afreq()

        # determine whether the allele frequency is between 0.0 and 1.0
        out = (afreq > 0.0) & (afreq < 1.0)

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
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A ``numpy.ndarray`` of shape ``(p,)`` containing allele frequencies for
            the minor allele.
        """
        out = self.afreq(dtype)     # get allele frequencies
        mask = out > 0.5            # create mask of allele frequencies > 0.5
        out[mask] = 1.0 - out[mask] # take 1 - allele frequency
        return out

    def meh(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> Real:
        """
        Mean expected heterozygosity across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : Real
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
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A ``numpy.ndarray`` array of shape ``(g,p)`` containing allele counts
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
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A ``numpy.ndarray`` array of shape ``(g,p)`` containing haplotype counts
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
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            overwrite: bool = True,
        ) -> None:
        """
        Write GenotypeMatrix to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name to which to write. File is closed after writing.
            If ``h5py.File``, an opened HDF5 file to which to write. File is not closed after writing.

        groupname : str, None
            If ``str``, an HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is written to the base HDF5 group.

        overwrite : bool
            Whether to overwrite values in an HDF5 file if a field already exists.
        """
        ########################################################
        ############ process ``filename`` argument #############

        # HDF5 file object
        h5file = None

        # if we have a string or Path, open HDF5 file in append (``r+``) mode
        if isinstance(filename, (str,Path)):
            h5file = h5py.File(filename, "a")

        # elif we have an h5py.File, make sure mode is writable, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_writable(filename)
            h5file = filename
        
        # else raise TypeError
        else:
            raise TypeError(
                "``filename`` must be of type ``str``, ``Path``, or ``h5py.File`` but received type ``{0}``".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        #### write data to HDF5 file and (optionally) close ####

        # data dictionary
        data = {
            "mat"               : self.mat,
            "taxa"              : self.taxa,
            "taxa_grp"          : self.taxa_grp,
            "vrnt_chrgrp"       : self.vrnt_chrgrp,
            "vrnt_phypos"       : self.vrnt_phypos,
            "vrnt_name"         : self.vrnt_name,
            "vrnt_genpos"       : self.vrnt_genpos,
            "vrnt_xoprob"       : self.vrnt_xoprob,
            "vrnt_hapgrp"       : self.vrnt_hapgrp,
            "vrnt_hapalt"       : self.vrnt_hapalt,
            "vrnt_hapref"       : self.vrnt_hapref,
            "vrnt_mask"         : self.vrnt_mask,
            "ploidy"            : self.ploidy,
            # metadata
            "taxa_grp_name"     : self.taxa_grp_name,
            "taxa_grp_stix"     : self.taxa_grp_stix,
            "taxa_grp_spix"     : self.taxa_grp_spix,
            "taxa_grp_len"      : self.taxa_grp_len,
            "vrnt_chrgrp_name"  : self.vrnt_chrgrp_name,
            "vrnt_chrgrp_stix"  : self.vrnt_chrgrp_stix,
            "vrnt_chrgrp_spix"  : self.vrnt_chrgrp_spix,
            "vrnt_chrgrp_len"   : self.vrnt_chrgrp_len,
        }

        # save data
        h5py_File_write_dict(h5file, groupname, data, overwrite)

        # close the file, only if the provided filename was a string or Path and not a h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None
        ) -> 'DenseGenotypeMatrix':
        """
        Read GenotypeMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which GenotypeMatrix data is stored.
            If ``None``, GenotypeMatrix is read from base HDF5 group.

        Returns
        -------
        gmat : GenotypeMatrix
            A genotype matrix read from file.
        """
        ########################################################
        ############ process ``filename`` argument #############

        # HDF5 file object
        h5file = None

        # if we have a string or Path, open HDF5 file in read (``r``) mode
        if isinstance(filename, (str,Path)):
            check_file_exists(filename)
            h5file = h5py.File(filename, "r")

        # elif we have an ``h5py.File``, make sure mode is in at least ``r`` mode, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_readable(filename)
            h5file = filename
        
        # else raise TypeError
        else:
            raise TypeError(
                "``filename`` must be of type ``str``, ``Path``, or ``h5py.File`` but received type ``{0}``".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # FIXME: errors if groupname == "" or "/"
            # if the group does not exist in the file, close and raise error
            check_h5py_File_has_group(h5file, groupname)

            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        ######## check that we have all required fields ########

        # all required arguments
        required_fields = ["mat"]

        # for each required field, check if the field exists in the HDF5 file.
        for field in required_fields:
            check_h5py_File_has_group(h5file, groupname + field)
        
        ########################################################
        ### read data from HDF5 file and (optionally) close ####
        
        # output dictionary
        data = {
            "mat"               : None,
            "taxa"              : None,
            "taxa_grp"          : None,
            "vrnt_chrgrp"       : None,
            "vrnt_phypos"       : None,
            "vrnt_name"         : None,
            "vrnt_genpos"       : None,
            "vrnt_xoprob"       : None,
            "vrnt_hapgrp"       : None,
            "vrnt_hapalt"       : None,
            "vrnt_hapref"       : None,
            "vrnt_mask"         : None,
            "ploidy"            : None,
            # metadata
            "taxa_grp_name"     : None,
            "taxa_grp_stix"     : None,
            "taxa_grp_spix"     : None,
            "taxa_grp_len"      : None,
            "vrnt_chrgrp_name"  : None,
            "vrnt_chrgrp_stix"  : None,
            "vrnt_chrgrp_spix"  : None,
            "vrnt_chrgrp_len"   : None,
        }

        ##################################
        ### read mandatory data fields ###

        # read mat array (ndarray dtype = int8)
        data["mat"] = h5py_File_read_ndarray_int8(h5file, groupname + "mat")
        
        #################################
        ### read optional data fields ###

        # read taxa array (ndarray dtype = unicode / object)
        if groupname + "taxa" in h5file:
            data["taxa"] = h5py_File_read_ndarray_utf8(h5file, groupname + "taxa")

        # read taxa_grp array (ndarray dtype = any)
        if groupname + "taxa_grp" in h5file:
            data["taxa_grp"] = h5py_File_read_ndarray(h5file, groupname + "taxa_grp")
        
        # read vrnt_chrgrp array (ndarray dtype = any)
        if groupname + "vrnt_chrgrp" in h5file:
            data["vrnt_chrgrp"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_chrgrp")

        # read vrnt_phypos array (ndarray dtype = any)
        if groupname + "vrnt_phypos" in h5file:
            data["vrnt_phypos"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_phypos")

        # read vrnt_name array (ndarray dtype = unicode / object)
        if groupname + "vrnt_name" in h5file:
            data["vrnt_name"] = h5py_File_read_ndarray_utf8(h5file, groupname + "vrnt_name")

        # read vrnt_genpos array (ndarray dtype = any)
        if groupname + "vrnt_genpos" in h5file:
            data["vrnt_genpos"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_genpos")

        # read vrnt_xoprob array (ndarray dtype = any)
        if groupname + "vrnt_xoprob" in h5file:
            data["vrnt_xoprob"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_xoprob")

        # read vrnt_hapgrp array (ndarray dtype = any)
        if groupname + "vrnt_hapgrp" in h5file:
            data["vrnt_hapgrp"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_hapgrp")

        # read vrnt_hapalt array (ndarray dtype = unicode / object)
        if groupname + "vrnt_hapalt" in h5file:
            data["vrnt_hapalt"] = h5py_File_read_ndarray_utf8(h5file, groupname + "vrnt_hapalt")

        # read vrnt_hapref array (ndarray dtype = unicode / object)
        if groupname + "vrnt_hapref" in h5file:
            data["vrnt_hapref"] = h5py_File_read_ndarray_utf8(h5file, groupname + "vrnt_hapref")

        # read vrnt_mask array (ndarray dtype = any)
        if groupname + "vrnt_mask" in h5file:
            data["vrnt_mask"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_mask")

        # read ploidy (dtype = int)
        if groupname + "ploidy" in h5file:
            data["ploidy"] = h5py_File_read_int(h5file, groupname + "ploidy")

        #####################################
        ### read optional metadata fields ###

        # read taxa_grp_name array (ndarray dtype = any)
        if groupname + "taxa_grp_name" in h5file:
            data["taxa_grp_name"] = h5py_File_read_ndarray(h5file, groupname + "taxa_grp_name")

        # read taxa_grp_stix array (ndarray dtype = any)
        if groupname + "taxa_grp_stix" in h5file:
            data["taxa_grp_stix"] = h5py_File_read_ndarray(h5file, groupname + "taxa_grp_stix")

        # read taxa_grp_spix array (ndarray dtype = any)
        if groupname + "taxa_grp_spix" in h5file:
            data["taxa_grp_spix"] = h5py_File_read_ndarray(h5file, groupname + "taxa_grp_spix")

        # read taxa_grp_len array (ndarray dtype = any)
        if groupname + "taxa_grp_len" in h5file:
            data["taxa_grp_len"] = h5py_File_read_ndarray(h5file, groupname + "taxa_grp_len")

        # read vrnt_chrgrp_name array (ndarray dtype = any)
        if groupname + "vrnt_chrgrp_name" in h5file:
            data["vrnt_chrgrp_name"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_chrgrp_name")

        # read vrnt_chrgrp_stix array (ndarray dtype = any)
        if groupname + "vrnt_chrgrp_stix" in h5file:
            data["vrnt_chrgrp_stix"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_chrgrp_stix")

        # read vrnt_chrgrp_spix array (ndarray dtype = any)
        if groupname + "vrnt_chrgrp_spix" in h5file:
            data["vrnt_chrgrp_spix"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_chrgrp_spix")

        # read vrnt_chrgrp_len array (ndarray dtype = any)
        if groupname + "vrnt_chrgrp_len" in h5file:
            data["vrnt_chrgrp_len"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_chrgrp_len")

        ######################
        ### close the file ###

        # close the file, only if the provided fieldname was a string or Path an not an h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

        ########################################################
        ################### Object creation ####################
        
        # create object from read data
        out = cls(
            mat         = data["mat"],
            taxa        = data["taxa"],
            taxa_grp    = data["taxa_grp"],
            vrnt_chrgrp = data["vrnt_chrgrp"],
            vrnt_phypos = data["vrnt_phypos"],
            vrnt_name   = data["vrnt_name"],
            vrnt_genpos = data["vrnt_genpos"],
            vrnt_xoprob = data["vrnt_xoprob"],
            vrnt_hapgrp = data["vrnt_hapgrp"],
            vrnt_hapalt = data["vrnt_hapalt"],
            vrnt_hapref = data["vrnt_hapref"],
            vrnt_mask   = data["vrnt_mask"], 
            ploidy      = data["ploidy"],
        )

        # copy metadata
        out.taxa_grp_name    = data["taxa_grp_name"]
        out.taxa_grp_stix    = data["taxa_grp_stix"]
        out.taxa_grp_spix    = data["taxa_grp_spix"]
        out.taxa_grp_len     = data["taxa_grp_len"]
        out.vrnt_chrgrp_name = data["vrnt_chrgrp_name"]
        out.vrnt_chrgrp_stix = data["vrnt_chrgrp_stix"]
        out.vrnt_chrgrp_spix = data["vrnt_chrgrp_spix"]
        out.vrnt_chrgrp_len  = data["vrnt_chrgrp_len"]

        return out

    @classmethod
    def from_vcf(
            cls, 
            filename: str,
            auto_group_vrnt: bool = True
        ) -> 'DenseGenotypeMatrix':
        """
        Create a DenseGenotypeMatrix from a VCF file.

        Parameters
        ----------
        filename : str
            Path to VCF file.
        auto_group_vrnt : bool
            Whether to group variants in returned genotype matrix.

        Returns
        -------
        out : DenseGenotypeMatrix
            An unphased genotype matrix with associated metadata from VCF file.
        """
        # make VCF iterator
        vcf = cyvcf2.VCF(filename)

        # extract taxa names from vcf header
        taxa = numpy.array(vcf.samples, dtype = object)

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
        vrnt_chrgrp = numpy.int64(vrnt_chrgrp)              # convert to int64 array
        vrnt_phypos = numpy.int64(vrnt_phypos)              # convert to int64 array
        vrnt_name = numpy.array(vrnt_name, dtype = object)  # convert to object array

        out = cls(
            mat = mat,
            taxa = taxa,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            ploidy = ploidy
        )

        if auto_group_vrnt:
            out.group_vrnt()

        return out



################################## Utilities ###################################
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
        raise TypeError("'{0}' must be a DenseGenotypeMatrix.".format(vname))
