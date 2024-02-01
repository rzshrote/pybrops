from numbers import Real
from typing import Optional, Sequence, Union
from numpy import ndarray
from numpy.typing import ArrayLike, DTypeLike
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.PhasedMatrix import PhasedMatrix
from pybrops.core.mat.TaxaMatrix import TaxaMatrix
from pybrops.core.mat.VariantMatrix import VariantMatrix
from pybrops.popgen.gmap.GeneticMap import GeneticMap
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.HaplotypeMatrix import HaplotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.gmat.PhasedHaplotypeMatrix import PhasedHaplotypeMatrix


class DummyGenotypeMatrix(GenotypeMatrix):
    def __add__(self, value):
        return super().__add__(value)
    def __and__(self, value):
        return super().__and__(value)
    def __copy__(self):
        return super().__copy__()
    def __deepcopy__(self, memo):
        return super().__deepcopy__(memo)
    def __delitem__(self, key):
        return super().__delitem__(key)
    def __divmod__(self, value):
        return super().__divmod__(value)
    def __eq__(self, value):
        return super().__eq__(value)
    def __floordiv__(self, value):
        return super().__floordiv__(value)
    def __ge__(self, value):
        return super().__ge__(value)
    def __getitem__(self, key):
        return super().__getitem__(key)
    def __gt__(self, value):
        return super().__gt__(value)
    def __iadd__(self, value):
        return super().__iadd__(value)
    def __iand__(self, value):
        return super().__iand__(value)
    def __ifloordiv__(self, value):
        return super().__ifloordiv__(value)
    def __ilshift__(self, value):
        return super().__ilshift__(value)
    def __imatmul__(self, value):
        return super().__imatmul__(value)
    def __imod__(self, value):
        return super().__imod__(value)
    def __imul__(self, value):
        return super().__imul__(value)
    def __ior__(self, value):
        return super().__ior__(value)
    def __ipow__(self, value):
        return super().__ipow__(value)
    def __irshift__(self, value):
        return super().__irshift__(value)
    def __isub__(self, value):
        return super().__isub__(value)
    def __iter__(self):
        return super().__iter__()
    def __itruediv__(self, value):
        return super().__itruediv__(value)
    def __ixor__(self, value):
        return super().__ixor__(value)
    def __le__(self, value):
        return super().__le__(value)
    def __len__(self):
        return super().__len__()
    def __lshift__(self, value):
        return super().__lshift__(value)
    def __lt__(self, value):
        return super().__lt__(value)
    def __matmul__(self, value):
        return super().__matmul__(value)
    def __mod__(self, value):
        return super().__mod__(value)
    def __mul__(self, value):
        return super().__mul__(value)
    def __ne__(self, value):
        return super().__ne__(value)
    def __or__(self, value):
        return super().__or__(value)
    def __pow__(self, value):
        return super().__pow__(value)
    def __radd__(self, value):
        return super().__radd__(value)
    def __rand__(self, value):
        return super().__rand__(value)
    def __rdivmod__(self, value):
        return super().__rdivmod__(value)
    def __rfloordiv__(self, value):
        return super().__rfloordiv__(value)
    def __rlshift__(self, value):
        return super().__rlshift__(value)
    def __rmatmul__(self, value):
        return super().__rmatmul__(value)
    def __rmod__(self, value):
        return super().__rmod__(value)
    def __rmul__(self, value):
        return super().__rmul__(value)
    def __ror__(self, value):
        return super().__ror__(value)
    def __rrshift__(self, value):
        return super().__rrshift__(value)
    def __rshift__(self, value):
        return super().__rshift__(value)
    def __rsub__(self, value):
        return super().__rsub__(value)
    def __rtruediv__(self, value):
        return super().__rtruediv__(value)
    def __rxor__(self, value):
        return super().__rxor__(value)
    def __setitem__(self, key, value):
        return super().__setitem__(key, value)
    def __sub__(self, value):
        return super().__sub__(value)
    def __truediv__(self, value):
        return super().__truediv__(value)
    def __xor__(self, value):
        return super().__xor__(value)
    def acount(self, dtype: DTypeLike) -> ndarray:
        return super().acount(dtype)
    def adjoin(self, values: Matrix | ndarray, axis: int, **kwargs: dict) -> Matrix:
        return super().adjoin(values, axis, **kwargs)
    def adjoin_taxa(self, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> TaxaMatrix:
        return super().adjoin_taxa(values, taxa, taxa_grp, **kwargs)
    def adjoin_vrnt(self, values: Matrix | ndarray, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray, vrnt_name: ndarray, vrnt_genpos: ndarray, vrnt_xoprob: ndarray, vrnt_hapgrp: ndarray, vrnt_mask: ndarray, **kwargs: dict) -> VariantMatrix:
        return super().adjoin_vrnt(values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs)
    def afreq(self, dtype: DTypeLike) -> ndarray:
        return super().afreq(dtype)
    def apoly(self, dtype: DTypeLike) -> ndarray:
        return super().apoly(dtype)
    def append(self, values: Matrix | ndarray, axis: int, **kwargs: dict) -> None:
        return super().append(values, axis, **kwargs)
    def append_taxa(self, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> None:
        return super().append_taxa(values, taxa, taxa_grp, **kwargs)
    def append_vrnt(self, values: Matrix | ndarray, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray, vrnt_name: ndarray, vrnt_genpos: ndarray, vrnt_xoprob: ndarray, vrnt_hapgrp: ndarray, vrnt_mask: ndarray, **kwargs: dict) -> None:
        return super().append_vrnt(values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs)
    @classmethod
    def concat(cls, mats: Sequence, axis: int, **kwargs: dict) -> Matrix:
        return super().concat(mats, axis, **kwargs)
    @classmethod
    def concat_taxa(cls, mats: Sequence, **kwargs: dict) -> TaxaMatrix:
        return super().concat_taxa(mats, **kwargs)
    @classmethod
    def concat_vrnt(cls, mats: Sequence, **kwargs: dict) -> VariantMatrix:
        return super().concat_vrnt(mats, **kwargs)
    def copy(self) -> Matrix:
        return super().copy()
    def deepcopy(self, memo: dict) -> Matrix:
        return super().deepcopy(memo)
    def delete(self, obj: int | slice | Sequence, axis: int, **kwargs: dict) -> Matrix:
        return super().delete(obj, axis, **kwargs)
    def delete_taxa(self, obj: int | slice | Sequence, **kwargs: dict) -> TaxaMatrix:
        return super().delete_taxa(obj, **kwargs)
    def delete_vrnt(self, obj: int | slice | Sequence, **kwargs: dict) -> VariantMatrix:
        return super().delete_vrnt(obj, **kwargs)
    @classmethod
    def from_hdf5(cls, filename: str, groupname: str | None) -> HDF5InputOutput:
        return super().from_hdf5(filename, groupname)
    def group(self, axis: int, **kwargs: dict) -> None:
        return super().group(axis, **kwargs)
    def group_taxa(self, **kwargs: dict) -> None:
        return super().group_taxa(**kwargs)
    def group_vrnt(self, **kwargs: dict) -> None:
        return super().group_vrnt(**kwargs)
    def gtcount(self, dtype: DTypeLike) -> ndarray:
        return super().gtcount(dtype)
    def gtfreq(self, dtype: DTypeLike) -> ndarray:
        return super().gtfreq(dtype)
    def incorp(self, obj: int | slice | Sequence, values: Matrix | ndarray, axis: int, **kwargs) -> None:
        return super().incorp(obj, values, axis, **kwargs)
    def incorp_taxa(self, obj: int | slice | Sequence, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> None:
        return super().incorp_taxa(obj, values, taxa, taxa_grp, **kwargs)
    def incorp_vrnt(self, obj: int | slice | Sequence, values: Matrix | ndarray, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray, vrnt_name: ndarray, vrnt_genpos: ndarray, vrnt_xoprob: ndarray, vrnt_hapgrp: ndarray, vrnt_mask: ndarray, **kwargs: dict) -> None:
        return super().incorp_vrnt(obj, values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs)
    def insert(self, obj: int | slice | Sequence, values: ArrayLike, axis: int, **kwargs: dict) -> Matrix:
        return super().insert(obj, values, axis, **kwargs)
    def insert_taxa(self, obj: int | slice | Sequence, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> TaxaMatrix:
        return super().insert_taxa(obj, values, taxa, taxa_grp, **kwargs)
    def insert_vrnt(self, obj: int | slice | Sequence, values: ArrayLike, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray, vrnt_name: ndarray, vrnt_genpos: ndarray, vrnt_xoprob: ndarray, vrnt_hapgrp: ndarray, vrnt_mask: ndarray, **kwargs: dict) -> VariantMatrix:
        return super().insert_vrnt(obj, values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs)
    def interp_genpos(self, gmap: GeneticMap, **kwargs: dict) -> None:
        return super().interp_genpos(gmap, **kwargs)
    def interp_xoprob(self, gmap: GeneticMap, gmapfn: GeneticMapFunction, **kwargs: dict) -> None:
        return super().interp_xoprob(gmap, gmapfn, **kwargs)
    def is_grouped(self, axis: int, **kwargs: dict) -> bool:
        return super().is_grouped(axis, **kwargs)
    def is_grouped_taxa(self, **kwargs: dict) -> bool:
        return super().is_grouped_taxa(**kwargs)
    def is_grouped_vrnt(self, **kwargs: dict) -> bool:
        return super().is_grouped_vrnt(**kwargs)
    def lexsort(self, keys: tuple | ndarray, axis: int, **kwargs: dict) -> ndarray:
        return super().lexsort(keys, axis, **kwargs)
    def lexsort_taxa(self, keys: tuple | ndarray, **kwargs: dict) -> ndarray:
        return super().lexsort_taxa(keys, **kwargs)
    def lexsort_vrnt(self, keys: tuple | ndarray, **kwargs: dict) -> ndarray:
        return super().lexsort_vrnt(keys, **kwargs)
    def maf(self, dtype: DTypeLike) -> ndarray:
        return super().maf(dtype)
    @property
    def mat(self) -> object:
        return super().mat
    @mat.setter
    def mat(self, value: object) -> None:
        super().mat = value
    def mat_asformat(self, format: str) -> ndarray:
        return super().mat_asformat(format)
    @property
    def mat_format(self) -> object:
        return super().mat_format
    @mat_format.setter
    def mat_format(self, value: object) -> None:
        super().mat_format = value
    @property
    def mat_ndim(self) -> object:
        return super().mat_ndim
    @mat_ndim.setter
    def mat_ndim(self, value: object) -> None:
        super().mat_ndim = value
    @property
    def mat_shape(self) -> object:
        return super().mat_shape
    @mat_shape.setter
    def mat_shape(self, value: object) -> None:
        super().mat_shape = value
    def meh(self, dtype: DTypeLike) -> Real:
        return super().meh(dtype)
    @property
    def nphase(self) -> object:
        return super().nphase
    @nphase.setter
    def nphase(self, value: object) -> None:
        super().nphase = value
    @property
    def ntaxa(self) -> object:
        return super().ntaxa
    @ntaxa.setter
    def ntaxa(self, value: object) -> None:
        super().ntaxa = value
    @property
    def nvrnt(self) -> object:
        return super().nvrnt
    @property
    def ploidy(self) -> object:
        return super().ploidy
    @ploidy.setter
    def ploidy(self, value: object) -> None:
        super().ploidy = value
    def remove(self, obj: int | slice | Sequence, axis: int, **kwargs: dict) -> None:
        return super().remove(obj, axis, **kwargs)
    def remove_taxa(self, obj: int | slice | Sequence, **kwargs: dict) -> None:
        return super().remove_taxa(obj, **kwargs)
    def remove_vrnt(self, obj: int | slice | Sequence, **kwargs: dict) -> None:
        return super().remove_vrnt(obj, **kwargs)
    def reorder(self, indices: ndarray | Sequence, axis: int, **kwargs: dict) -> None:
        return super().reorder(indices, axis, **kwargs)
    def reorder_taxa(self, indices: ndarray | Sequence, **kwargs: dict) -> None:
        return super().reorder_taxa(indices, **kwargs)
    def reorder_vrnt(self, indices: ndarray | Sequence, **kwargs: dict) -> None:
        return super().reorder_vrnt(indices, **kwargs)
    def select(self, indices: ArrayLike, axis: int, **kwargs: dict) -> Matrix:
        return super().select(indices, axis, **kwargs)
    def select_taxa(self, indices: ArrayLike, **kwargs: dict) -> TaxaMatrix:
        return super().select_taxa(indices, **kwargs)
    def select_vrnt(self, indices: ArrayLike, **kwargs: dict) -> VariantMatrix:
        return super().select_vrnt(indices, **kwargs)
    def sort(self, keys: tuple | ndarray, axis: int, **kwargs: dict) -> None:
        return super().sort(keys, axis, **kwargs)
    def sort_taxa(self, keys: tuple | ndarray, **kwargs: dict) -> None:
        return super().sort_taxa(keys, **kwargs)
    def sort_vrnt(self, keys: tuple | ndarray, **kwargs: dict) -> None:
        return super().sort_vrnt(keys, **kwargs)
    def tacount(self, dtype: DTypeLike) -> ndarray:
        return super().tacount(dtype)
    def tafreq(self, dtype: DTypeLike) -> ndarray:
        return super().tafreq(dtype)
    @property
    def taxa(self) -> object:
        return super().taxa
    @taxa.setter
    def taxa(self, value: object) -> None:
        super().taxa = value
    @property
    def taxa_axis(self) -> object:
        return super().taxa_axis
    @taxa_axis.setter
    def taxa_axis(self, value: object) -> None:
        super().taxa_axis = value
    @property
    def taxa_grp(self) -> object:
        return super().taxa_grp
    @taxa_grp.setter
    def taxa_grp(self, value: object) -> None:
        super().taxa_grp = value
    @property
    def taxa_grp_len(self) -> object:
        return super().taxa_grp_len
    @taxa_grp_len.setter
    def taxa_grp_len(self, value: object) -> None:
        super().taxa_grp_len = value
    @property
    def taxa_grp_name(self) -> object:
        return super().taxa_grp_name
    @taxa_grp_name.setter
    def taxa_grp_name(self, value: object) -> None:
        super().taxa_grp_name = value
    @property
    def taxa_grp_stix(self) -> object:
        return super().taxa_grp_stix
    @taxa_grp_stix.setter
    def taxa_grp_stix(self, value: object) -> None:
        super().taxa_grp_stix = value
    @property
    def taxa_grp_spix(self) -> object:
        return super().taxa_grp_spix
    @taxa_grp_spix.setter
    def taxa_grp_spix(self, value: object) -> None:
        super().taxa_grp_spix = value
    def to_hdf5(self, filename: str, groupname: str | None) -> None:
        return super().to_hdf5(filename, groupname)
    @property
    def vrnt_axis(self) -> object:
        return super().vrnt_axis
    @vrnt_axis.setter
    def vrnt_axis(self, value: object) -> None:
        super().vrnt_axis = value
    @property
    def vrnt_chrgrp(self) -> object:
        return super().vrnt_chrgrp
    @vrnt_chrgrp.setter
    def vrnt_chrgrp(self, value: object) -> None:
        super().vrnt_chrgrp = value
    @property
    def vrnt_chrgrp_len(self) -> object:
        return super().vrnt_chrgrp_len
    @vrnt_chrgrp_len.setter
    def vrnt_chrgrp_len(self, value: object) -> None:
        super().vrnt_chrgrp_len = value
    @property
    def vrnt_chrgrp_name(self) -> object:
        return super().vrnt_chrgrp_name
    @vrnt_chrgrp_name.setter
    def vrnt_chrgrp_name(self, value: object) -> None:
        super().vrnt_chrgrp_name = value
    @property
    def vrnt_chrgrp_stix(self) -> object:
        return super().vrnt_chrgrp_stix
    @vrnt_chrgrp_stix.setter
    def vrnt_chrgrp_stix(self, value: object) -> None:
        super().vrnt_chrgrp_stix = value
    @property
    def vrnt_chrgrp_spix(self) -> object:
        return super().vrnt_chrgrp_spix
    @vrnt_chrgrp_spix.setter
    def vrnt_chrgrp_spix(self, value: object) -> None:
        super().vrnt_chrgrp_spix = value
    @property
    def vrnt_genpos(self) -> object:
        return super().vrnt_genpos
    @vrnt_genpos.setter
    def vrnt_genpos(self, value: object) -> None:
        super().vrnt_genpos = value
    @property
    def vrnt_hapalt(self) -> object:
        return super().vrnt_hapalt
    @vrnt_hapalt.setter
    def vrnt_hapalt(self, value: object) -> None:
        super().vrnt_hapalt = value
    @property
    def vrnt_hapref(self) -> object:
        return super().vrnt_hapref
    @vrnt_hapref.setter
    def vrnt_hapref(self, value: object) -> None:
        super().vrnt_hapref = value
    @property
    def vrnt_hapgrp(self) -> object:
        return super().vrnt_hapgrp
    @vrnt_hapgrp.setter
    def vrnt_hapgrp(self, value: object) -> None:
        super().vrnt_hapgrp = value
    @property
    def vrnt_mask(self) -> object:
        return super().vrnt_mask
    @vrnt_mask.setter
    def vrnt_mask(self, value: object) -> None:
        super().vrnt_mask = value
    @property
    def vrnt_name(self) -> object:
        return super().vrnt_name
    @vrnt_name.setter
    def vrnt_name(self, value: object) -> None:
        super().vrnt_name = value
    @property
    def vrnt_phypos(self) -> object:
        return super().vrnt_phypos
    @vrnt_phypos.setter
    def vrnt_phypos(self, value: object) -> None:
        super().vrnt_phypos = value
    @property
    def vrnt_xoprob(self) -> object:
        return super().vrnt_xoprob
    @vrnt_xoprob.setter
    def vrnt_xoprob(self, value: object) -> None:
        super().vrnt_xoprob = value
    def ungroup(self, axis: int, **kwargs: dict) -> None:
        return super().ungroup(axis, **kwargs)
    def ungroup_taxa(self, **kwargs: dict) -> None:
        return super().ungroup_taxa(**kwargs)
    def ungroup_vrnt(self, **kwargs: dict) -> None:
        return super().ungroup_vrnt(**kwargs)

class DummyHaplotypeMatrix(HaplotypeMatrix):
    def __add__(self, value):
        return super().__add__(value)
    def __and__(self, value):
        return super().__and__(value)
    def __copy__(self):
        return super().__copy__()
    def __deepcopy__(self, memo):
        return super().__deepcopy__(memo)
    def __delitem__(self, key):
        return super().__delitem__(key)
    def __divmod__(self, value):
        return super().__divmod__(value)
    def __eq__(self, value):
        return super().__eq__(value)
    def __floordiv__(self, value):
        return super().__floordiv__(value)
    def __ge__(self, value):
        return super().__ge__(value)
    def __getitem__(self, key):
        return super().__getitem__(key)
    def __gt__(self, value):
        return super().__gt__(value)
    def __iadd__(self, value):
        return super().__iadd__(value)
    def __iand__(self, value):
        return super().__iand__(value)
    def __ifloordiv__(self, value):
        return super().__ifloordiv__(value)
    def __ilshift__(self, value):
        return super().__ilshift__(value)
    def __imatmul__(self, value):
        return super().__imatmul__(value)
    def __imod__(self, value):
        return super().__imod__(value)
    def __imul__(self, value):
        return super().__imul__(value)
    def __ior__(self, value):
        return super().__ior__(value)
    def __ipow__(self, value):
        return super().__ipow__(value)
    def __irshift__(self, value):
        return super().__irshift__(value)
    def __isub__(self, value):
        return super().__isub__(value)
    def __iter__(self):
        return super().__iter__()
    def __itruediv__(self, value):
        return super().__itruediv__(value)
    def __ixor__(self, value):
        return super().__ixor__(value)
    def __le__(self, value):
        return super().__le__(value)
    def __len__(self):
        return super().__len__()
    def __lshift__(self, value):
        return super().__lshift__(value)
    def __lt__(self, value):
        return super().__lt__(value)
    def __matmul__(self, value):
        return super().__matmul__(value)
    def __mod__(self, value):
        return super().__mod__(value)
    def __mul__(self, value):
        return super().__mul__(value)
    def __ne__(self, value):
        return super().__ne__(value)
    def __or__(self, value):
        return super().__or__(value)
    def __pow__(self, value):
        return super().__pow__(value)
    def __radd__(self, value):
        return super().__radd__(value)
    def __rand__(self, value):
        return super().__rand__(value)
    def __rdivmod__(self, value):
        return super().__rdivmod__(value)
    def __rfloordiv__(self, value):
        return super().__rfloordiv__(value)
    def __rlshift__(self, value):
        return super().__rlshift__(value)
    def __rmatmul__(self, value):
        return super().__rmatmul__(value)
    def __rmod__(self, value):
        return super().__rmod__(value)
    def __rmul__(self, value):
        return super().__rmul__(value)
    def __ror__(self, value):
        return super().__ror__(value)
    def __rrshift__(self, value):
        return super().__rrshift__(value)
    def __rshift__(self, value):
        return super().__rshift__(value)
    def __rsub__(self, value):
        return super().__rsub__(value)
    def __rtruediv__(self, value):
        return super().__rtruediv__(value)
    def __rxor__(self, value):
        return super().__rxor__(value)
    def __setitem__(self, key, value):
        return super().__setitem__(key, value)
    def __sub__(self, value):
        return super().__sub__(value)
    def __truediv__(self, value):
        return super().__truediv__(value)
    def __xor__(self, value):
        return super().__xor__(value)
    def adjoin(self, values: Matrix | ndarray, axis: int, **kwargs: dict) -> Matrix:
        return super().adjoin(values, axis, **kwargs)
    def adjoin_taxa(self, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> TaxaMatrix:
        return super().adjoin_taxa(values, taxa, taxa_grp, **kwargs)
    def adjoin_vrnt(self, values: Matrix | ndarray, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray, vrnt_name: ndarray, vrnt_genpos: ndarray, vrnt_xoprob: ndarray, vrnt_hapgrp: ndarray, vrnt_mask: ndarray, **kwargs: dict) -> VariantMatrix:
        return super().adjoin_vrnt(values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs)
    def append(self, values: Matrix | ndarray, axis: int, **kwargs: dict) -> None:
        return super().append(values, axis, **kwargs)
    def append_taxa(self, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> None:
        return super().append_taxa(values, taxa, taxa_grp, **kwargs)
    def append_vrnt(self, values: Matrix | ndarray, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray, vrnt_name: ndarray, vrnt_genpos: ndarray, vrnt_xoprob: ndarray, vrnt_hapgrp: ndarray, vrnt_mask: ndarray, **kwargs: dict) -> None:
        return super().append_vrnt(values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs)
    @classmethod
    def concat(cls, mats: Sequence, axis: int, **kwargs: dict) -> Matrix:
        return super().concat(mats, axis, **kwargs)
    @classmethod
    def concat_taxa(cls, mats: Sequence, **kwargs: dict) -> TaxaMatrix:
        return super().concat_taxa(mats, **kwargs)
    @classmethod
    def concat_vrnt(cls, mats: Sequence, **kwargs: dict) -> VariantMatrix:
        return super().concat_vrnt(mats, **kwargs)
    def copy(self) -> Matrix:
        return super().copy()
    def deepcopy(self, memo: dict) -> Matrix:
        return super().deepcopy(memo)
    def delete(self, obj: int | slice | Sequence, axis: int, **kwargs: dict) -> Matrix:
        return super().delete(obj, axis, **kwargs)
    def delete_taxa(self, obj: int | slice | Sequence, **kwargs: dict) -> TaxaMatrix:
        return super().delete_taxa(obj, **kwargs)
    def delete_vrnt(self, obj: int | slice | Sequence, **kwargs: dict) -> VariantMatrix:
        return super().delete_vrnt(obj, **kwargs)
    @classmethod
    def from_hdf5(cls, filename: str, groupname: str | None) -> HDF5InputOutput:
        return super().from_hdf5(filename, groupname)
    def group(self, axis: int, **kwargs: dict) -> None:
        return super().group(axis, **kwargs)
    def group_taxa(self, **kwargs: dict) -> None:
        return super().group_taxa(**kwargs)
    def group_vrnt(self, **kwargs: dict) -> None:
        return super().group_vrnt(**kwargs)
    def gtcount(self, dtype: DTypeLike) -> ndarray:
        return super().gtcount(dtype)
    def gtfreq(self, dtype: DTypeLike) -> ndarray:
        return super().gtfreq(dtype)
    def incorp(self, obj: int | slice | Sequence, values: Matrix | ndarray, axis: int, **kwargs) -> None:
        return super().incorp(obj, values, axis, **kwargs)
    def incorp_taxa(self, obj: int | slice | Sequence, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> None:
        return super().incorp_taxa(obj, values, taxa, taxa_grp, **kwargs)
    def incorp_vrnt(self, obj: int | slice | Sequence, values: Matrix | ndarray, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray, vrnt_name: ndarray, vrnt_genpos: ndarray, vrnt_xoprob: ndarray, vrnt_hapgrp: ndarray, vrnt_mask: ndarray, **kwargs: dict) -> None:
        return super().incorp_vrnt(obj, values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs)
    def insert(self, obj: int | slice | Sequence, values: ArrayLike, axis: int, **kwargs: dict) -> Matrix:
        return super().insert(obj, values, axis, **kwargs)
    def insert_taxa(self, obj: int | slice | Sequence, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> TaxaMatrix:
        return super().insert_taxa(obj, values, taxa, taxa_grp, **kwargs)
    def insert_vrnt(self, obj: int | slice | Sequence, values: ArrayLike, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray, vrnt_name: ndarray, vrnt_genpos: ndarray, vrnt_xoprob: ndarray, vrnt_hapgrp: ndarray, vrnt_mask: ndarray, **kwargs: dict) -> VariantMatrix:
        return super().insert_vrnt(obj, values, vrnt_chrgrp, vrnt_phypos, vrnt_name, vrnt_genpos, vrnt_xoprob, vrnt_hapgrp, vrnt_mask, **kwargs)
    def is_grouped(self, axis: int, **kwargs: dict) -> bool:
        return super().is_grouped(axis, **kwargs)
    def is_grouped_taxa(self, **kwargs: dict) -> bool:
        return super().is_grouped_taxa(**kwargs)
    def is_grouped_vrnt(self, **kwargs: dict) -> bool:
        return super().is_grouped_vrnt(**kwargs)
    def lexsort(self, keys: tuple | ndarray, axis: int, **kwargs: dict) -> ndarray:
        return super().lexsort(keys, axis, **kwargs)
    def lexsort_taxa(self, keys: tuple | ndarray, **kwargs: dict) -> ndarray:
        return super().lexsort_taxa(keys, **kwargs)
    def lexsort_vrnt(self, keys: tuple | ndarray, **kwargs: dict) -> ndarray:
        return super().lexsort_vrnt(keys, **kwargs)
    @property
    def mat(self) -> object:
        return super().mat
    @mat.setter
    def mat(self, value: object) -> None:
        super().mat = value
    @property
    def mat_format(self) -> object:
        return super().mat_format
    @mat_format.setter
    def mat_format(self, value: object) -> None:
        super().mat_format = value
    @property
    def mat_ndim(self) -> object:
        return super().mat_ndim
    @mat_ndim.setter
    def mat_ndim(self, value: object) -> None:
        super().mat_ndim = value
    @property
    def mat_shape(self) -> object:
        return super().mat_shape
    @mat_shape.setter
    def mat_shape(self, value: object) -> None:
        super().mat_shape = value
    def meh(self, dtype: DTypeLike) -> Real:
        return super().meh(dtype)
    @property
    def nphase(self) -> object:
        return super().nphase
    @nphase.setter
    def nphase(self, value: object) -> None:
        super().nphase = value
    @property
    def ntaxa(self) -> object:
        return super().ntaxa
    @ntaxa.setter
    def ntaxa(self, value: object) -> None:
        super().ntaxa = value
    @property
    def nvrnt(self) -> object:
        return super().nvrnt
    @property
    def ploidy(self) -> object:
        return super().ploidy
    @ploidy.setter
    def ploidy(self, value: object) -> None:
        super().ploidy = value
    def remove(self, obj: int | slice | Sequence, axis: int, **kwargs: dict) -> None:
        return super().remove(obj, axis, **kwargs)
    def remove_taxa(self, obj: int | slice | Sequence, **kwargs: dict) -> None:
        return super().remove_taxa(obj, **kwargs)
    def remove_vrnt(self, obj: int | slice | Sequence, **kwargs: dict) -> None:
        return super().remove_vrnt(obj, **kwargs)
    def reorder(self, indices: ndarray | Sequence, axis: int, **kwargs: dict) -> None:
        return super().reorder(indices, axis, **kwargs)
    def reorder_taxa(self, indices: ndarray | Sequence, **kwargs: dict) -> None:
        return super().reorder_taxa(indices, **kwargs)
    def reorder_vrnt(self, indices: ndarray | Sequence, **kwargs: dict) -> None:
        return super().reorder_vrnt(indices, **kwargs)
    def select(self, indices: ArrayLike, axis: int, **kwargs: dict) -> Matrix:
        return super().select(indices, axis, **kwargs)
    def select_taxa(self, indices: ArrayLike, **kwargs: dict) -> TaxaMatrix:
        return super().select_taxa(indices, **kwargs)
    def select_vrnt(self, indices: ArrayLike, **kwargs: dict) -> VariantMatrix:
        return super().select_vrnt(indices, **kwargs)
    def sort(self, keys: tuple | ndarray, axis: int, **kwargs: dict) -> None:
        return super().sort(keys, axis, **kwargs)
    def sort_taxa(self, keys: tuple | ndarray, **kwargs: dict) -> None:
        return super().sort_taxa(keys, **kwargs)
    def sort_vrnt(self, keys: tuple | ndarray, **kwargs: dict) -> None:
        return super().sort_vrnt(keys, **kwargs)
    @property
    def taxa(self) -> object:
        return super().taxa
    @taxa.setter
    def taxa(self, value: object) -> None:
        super().taxa = value
    @property
    def taxa_axis(self) -> object:
        return super().taxa_axis
    @taxa_axis.setter
    def taxa_axis(self, value: object) -> None:
        super().taxa_axis = value
    @property
    def taxa_grp(self) -> object:
        return super().taxa_grp
    @taxa_grp.setter
    def taxa_grp(self, value: object) -> None:
        super().taxa_grp = value
    @property
    def taxa_grp_len(self) -> object:
        return super().taxa_grp_len
    @taxa_grp_len.setter
    def taxa_grp_len(self, value: object) -> None:
        super().taxa_grp_len = value
    @property
    def taxa_grp_name(self) -> object:
        return super().taxa_grp_name
    @taxa_grp_name.setter
    def taxa_grp_name(self, value: object) -> None:
        super().taxa_grp_name = value
    @property
    def taxa_grp_stix(self) -> object:
        return super().taxa_grp_stix
    @taxa_grp_stix.setter
    def taxa_grp_stix(self, value: object) -> None:
        super().taxa_grp_stix = value
    @property
    def taxa_grp_spix(self) -> object:
        return super().taxa_grp_spix
    @taxa_grp_spix.setter
    def taxa_grp_spix(self, value: object) -> None:
        super().taxa_grp_spix = value
    def to_hdf5(self, filename: str, groupname: str | None) -> None:
        return super().to_hdf5(filename, groupname)
    @property
    def vrnt_axis(self) -> object:
        return super().vrnt_axis
    @vrnt_axis.setter
    def vrnt_axis(self, value: object) -> None:
        super().vrnt_axis = value
    @property
    def vrnt_chrgrp(self) -> object:
        return super().vrnt_chrgrp
    @vrnt_chrgrp.setter
    def vrnt_chrgrp(self, value: object) -> None:
        super().vrnt_chrgrp = value
    @property
    def vrnt_chrgrp_len(self) -> object:
        return super().vrnt_chrgrp_len
    @vrnt_chrgrp_len.setter
    def vrnt_chrgrp_len(self, value: object) -> None:
        super().vrnt_chrgrp_len = value
    @property
    def vrnt_chrgrp_name(self) -> object:
        return super().vrnt_chrgrp_name
    @vrnt_chrgrp_name.setter
    def vrnt_chrgrp_name(self, value: object) -> None:
        super().vrnt_chrgrp_name = value
    @property
    def vrnt_chrgrp_stix(self) -> object:
        return super().vrnt_chrgrp_stix
    @vrnt_chrgrp_stix.setter
    def vrnt_chrgrp_stix(self, value: object) -> None:
        super().vrnt_chrgrp_stix = value
    @property
    def vrnt_chrgrp_spix(self) -> object:
        return super().vrnt_chrgrp_spix
    @vrnt_chrgrp_spix.setter
    def vrnt_chrgrp_spix(self, value: object) -> None:
        super().vrnt_chrgrp_spix = value
    @property
    def vrnt_genpos(self) -> object:
        return super().vrnt_genpos
    @vrnt_genpos.setter
    def vrnt_genpos(self, value: object) -> None:
        super().vrnt_genpos = value
    @property
    def vrnt_hapalt(self) -> object:
        return super().vrnt_hapalt
    @vrnt_hapalt.setter
    def vrnt_hapalt(self, value: object) -> None:
        super().vrnt_hapalt = value
    @property
    def vrnt_hapref(self) -> object:
        return super().vrnt_hapref
    @vrnt_hapref.setter
    def vrnt_hapref(self, value: object) -> None:
        super().vrnt_hapref = value
    @property
    def vrnt_hapgrp(self) -> object:
        return super().vrnt_hapgrp
    @vrnt_hapgrp.setter
    def vrnt_hapgrp(self, value: object) -> None:
        super().vrnt_hapgrp = value
    @property
    def vrnt_mask(self) -> object:
        return super().vrnt_mask
    @vrnt_mask.setter
    def vrnt_mask(self, value: object) -> None:
        super().vrnt_mask = value
    @property
    def vrnt_name(self) -> object:
        return super().vrnt_name
    @vrnt_name.setter
    def vrnt_name(self, value: object) -> None:
        super().vrnt_name = value
    @property
    def vrnt_phypos(self) -> object:
        return super().vrnt_phypos
    @vrnt_phypos.setter
    def vrnt_phypos(self, value: object) -> None:
        super().vrnt_phypos = value
    @property
    def vrnt_xoprob(self) -> object:
        return super().vrnt_xoprob
    @vrnt_xoprob.setter
    def vrnt_xoprob(self, value: object) -> None:
        super().vrnt_xoprob = value
    @property
    def nphase(self) -> object:
        return super().nphase
    @nphase.setter
    def nphase(self, value: object) -> None:
        super().nphase = value
    def hcount(self, dtype: DTypeLike) -> ndarray:
        return super().hcount(dtype)
    def hfreq(self, dtype: DTypeLike) -> ndarray:
        return super().hfreq(dtype)
    def mhf(self, dtype: DTypeLike) -> ndarray:
        return super().mhf(dtype)
    def thcount(self, dtype: DTypeLike) -> ndarray:
        return super().thcount(dtype)
    def thfreq(self, dtype: DTypeLike) -> ndarray:
        return super().thfreq(dtype)
    def ungroup(self, axis: int, **kwargs: dict) -> None:
        return super().ungroup(axis, **kwargs)
    def ungroup_taxa(self, **kwargs: dict) -> None:
        return super().ungroup_taxa(**kwargs)
    def ungroup_vrnt(self, **kwargs: dict) -> None:
        return super().ungroup_vrnt(**kwargs)

class DummyPhasedGenotypeMatrix(DummyGenotypeMatrix,PhasedGenotypeMatrix):
    def adjoin_phase(self, values: Matrix | ndarray, **kwargs: dict) -> PhasedMatrix:
        return super().adjoin_phase(values, **kwargs)
    def append_phase(self, values: Matrix | ndarray, **kwargs: dict) -> None:
        return super().append_phase(values, **kwargs)
    @classmethod
    def concat_phase(cls, mats: Sequence, **kwargs: dict) -> PhasedMatrix:
        return super().concat_phase(mats, **kwargs)
    def delete_phase(self, obj: int | slice | Sequence, **kwargs: dict) -> PhasedMatrix:
        return super().delete_phase(obj, **kwargs)
    def incorp_phase(self, obj: int | slice | Sequence, values: Matrix | ndarray, **kwargs: dict) -> None:
        return super().incorp_phase(obj, values, **kwargs)
    def insert_phase(self, obj: int | slice | Sequence, values: Matrix | ndarray, **kwargs: dict) -> PhasedMatrix:
        return super().insert_phase(obj, values, **kwargs)
    @property
    def phase_axis(self) -> object:
        return self._phase_axis
    @phase_axis.setter
    def phase_axis(self, value: object) -> None:
        self._phase_axis = value
    def remove_phase(self, obj: int | slice | Sequence, **kwargs: dict) -> None:
        return super().remove_phase(obj, **kwargs)
    def select_phase(self, indices: ArrayLike, **kwargs: dict) -> PhasedMatrix:
        return super().select_phase(indices, **kwargs)
    
class DummyPhasedHaplotypeMatrix(DummyHaplotypeMatrix,PhasedHaplotypeMatrix):
    def adjoin_phase(self, values: Matrix | ndarray, **kwargs: dict) -> PhasedMatrix:
        return super().adjoin_phase(values, **kwargs)
    def append_phase(self, values: Matrix | ndarray, **kwargs: dict) -> None:
        return super().append_phase(values, **kwargs)
    @classmethod
    def concat_phase(cls, mats: Sequence, **kwargs: dict) -> PhasedMatrix:
        return super().concat_phase(mats, **kwargs)
    def delete_phase(self, obj: int | slice | Sequence, **kwargs: dict) -> PhasedMatrix:
        return super().delete_phase(obj, **kwargs)
    def incorp_phase(self, obj: int | slice | Sequence, values: Matrix | ndarray, **kwargs: dict) -> None:
        return super().incorp_phase(obj, values, **kwargs)
    def insert_phase(self, obj: int | slice | Sequence, values: Matrix | ndarray, **kwargs: dict) -> PhasedMatrix:
        return super().insert_phase(obj, values, **kwargs)
    @property
    def phase_axis(self) -> object:
        return self._phase_axis
    @phase_axis.setter
    def phase_axis(self, value: object) -> None:
        self._phase_axis = value
    def remove_phase(self, obj: int | slice | Sequence, **kwargs: dict) -> None:
        return super().remove_phase(obj, **kwargs)
    def select_phase(self, indices: ArrayLike, **kwargs: dict) -> PhasedMatrix:
        return super().select_phase(indices, **kwargs)
    