from numbers import Integral
from typing import Sequence

import numpy
from numpy.typing import ArrayLike
import pandas
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.TaxaMatrix import TaxaMatrix
from pybrops.core.mat.TraitMatrix import TraitMatrix
from pybrops.model.embvmat.ExpectedMaximumBreedingValueMatrix import ExpectedMaximumBreedingValueMatrix
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class DummyExpectedMaximumBreedingValueMatrix(ExpectedMaximumBreedingValueMatrix):
    def __init__(self) -> None:
        pass
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
    def adjoin(self, values: Matrix | numpy.ndarray, axis: int, **kwargs: dict) -> Matrix:
        return super().adjoin(values, axis, **kwargs)
    def adjoin_taxa(self, values: Matrix | numpy.ndarray, taxa: numpy.ndarray, taxa_grp: numpy.ndarray, **kwargs: dict) -> TaxaMatrix:
        return super().adjoin_taxa(values, taxa, taxa_grp, **kwargs)
    def adjoin_trait(self, values: Matrix | numpy.ndarray, trait: numpy.ndarray, **kwargs: dict) -> TraitMatrix:
        return super().adjoin_trait(values, trait, **kwargs)
    def append(self, values: Matrix | numpy.ndarray, axis: int, **kwargs: dict) -> None:
        return super().append(values, axis, **kwargs)
    def append_taxa(self, values: Matrix | numpy.ndarray, taxa: numpy.ndarray, taxa_grp: numpy.ndarray, **kwargs: dict) -> None:
        return super().append_taxa(values, taxa, taxa_grp, **kwargs)
    def append_trait(self, values: Matrix | numpy.ndarray, trait: numpy.ndarray, **kwargs: dict) -> None:
        return super().append_trait(values, trait, **kwargs)
    @classmethod
    def concat(cls, mats: Sequence, axis: int, **kwargs: dict) -> Matrix:
        return super().concat(mats, axis, **kwargs)
    @classmethod
    def concat_taxa(cls, mats: Sequence, **kwargs: dict) -> TaxaMatrix:
        return super().concat_taxa(mats, **kwargs)
    @classmethod
    def concat_trait(cls, mats: Sequence, **kwargs: dict) -> TraitMatrix:
        return super().concat_trait(mats, **kwargs)
    def copy(self) -> Matrix:
        return super().copy()
    def deepcopy(self, memo: dict) -> Matrix:
        return super().deepcopy(memo)
    def delete(self, obj: int | slice | Sequence, axis: int, **kwargs: dict) -> Matrix:
        return super().delete(obj, axis, **kwargs)
    def delete_taxa(self, obj: int | slice | Sequence, **kwargs: dict) -> TaxaMatrix:
        return super().delete_taxa(obj, **kwargs)
    def delete_trait(self, obj: int | slice | Sequence, **kwargs: dict) -> TraitMatrix:
        return super().delete_trait(obj, **kwargs)
    @classmethod
    def from_hdf5(cls, filename: str, groupname: str | None) -> HDF5InputOutput:
        return super().from_hdf5(filename, groupname)
    def group(self, axis: int, **kwargs: dict) -> None:
        return super().group(axis, **kwargs)
    def group_taxa(self, **kwargs: dict) -> None:
        return super().group_taxa(**kwargs)
    def incorp(self, obj: int | slice | Sequence, values: Matrix | numpy.ndarray, axis: int, **kwargs) -> None:
        return super().incorp(obj, values, axis, **kwargs)
    def incorp_taxa(self, obj: int | slice | Sequence, values: Matrix | numpy.ndarray, taxa: numpy.ndarray, taxa_grp: numpy.ndarray, **kwargs: dict) -> None:
        return super().incorp_taxa(obj, values, taxa, taxa_grp, **kwargs)
    def incorp_trait(self, obj: int | slice | Sequence, values: Matrix | numpy.ndarray, trait: numpy.ndarray, **kwargs: dict) -> None:
        return super().incorp_trait(obj, values, trait, **kwargs)
    def insert(self, obj: int | slice | Sequence, values: ArrayLike, axis: int, **kwargs: dict) -> Matrix:
        return super().insert(obj, values, axis, **kwargs)
    def insert_taxa(self, obj: int | slice | Sequence, values: Matrix | numpy.ndarray, taxa: numpy.ndarray, taxa_grp: numpy.ndarray, **kwargs: dict) -> TaxaMatrix:
        return super().insert_taxa(obj, values, taxa, taxa_grp, **kwargs)
    def insert_trait(self, obj: int | slice | Sequence, values: Matrix | numpy.ndarray, trait: numpy.ndarray, **kwargs) -> TraitMatrix:
        return super().insert_trait(obj, values, trait, **kwargs)
    def is_grouped(self, axis: int, **kwargs: dict) -> bool:
        return super().is_grouped(axis, **kwargs)
    def is_grouped_taxa(self, **kwargs: dict) -> bool:
        return super().is_grouped_taxa(**kwargs)
    def lexsort(self, keys: tuple | numpy.ndarray, axis: int, **kwargs: dict) -> numpy.ndarray:
        return super().lexsort(keys, axis, **kwargs)
    def lexsort_taxa(self, keys: tuple | numpy.ndarray, **kwargs: dict) -> numpy.ndarray:
        return super().lexsort_taxa(keys, **kwargs)
    def lexsort_trait(self, keys: tuple | numpy.ndarray, **kwargs: dict) -> numpy.ndarray:
        return super().lexsort_trait(keys, **kwargs)
    @property
    def mat(self) -> object:
        return super().mat
    @mat.setter
    def mat(self, value: object) -> None:
        super().mat = value
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
    @property
    def ntaxa(self) -> object:
        return super().ntaxa
    @ntaxa.setter
    def ntaxa(self, value: object) -> None:
        super().ntaxa = value
    @property
    def ntrait(self) -> object:
        return super().ntrait
    @ntrait.setter
    def ntrait(self, value: object) -> None:
        super().ntrait = value
    def remove(self, obj: int | slice | Sequence, axis: int, **kwargs: dict) -> None:
        return super().remove(obj, axis, **kwargs)
    def remove_taxa(self, obj: int | slice | Sequence, **kwargs: dict) -> None:
        return super().remove_taxa(obj, **kwargs)
    def remove_trait(self, obj: int | slice | Sequence, **kwargs: dict) -> None:
        return super().remove_trait(obj, **kwargs)
    def reorder(self, indices: numpy.ndarray | Sequence, axis: int, **kwargs: dict) -> None:
        return super().reorder(indices, axis, **kwargs)
    def reorder_taxa(self, indices: numpy.ndarray | Sequence, **kwargs: dict) -> None:
        return super().reorder_taxa(indices, **kwargs)
    def reorder_trait(self, indices: numpy.ndarray | Sequence, **kwargs: dict) -> None:
        return super().reorder_trait(indices, **kwargs)
    def select(self, indices: ArrayLike, axis: int, **kwargs: dict) -> Matrix:
        return super().select(indices, axis, **kwargs)
    def select_taxa(self, indices: ArrayLike, **kwargs: dict) -> TaxaMatrix:
        return super().select_taxa(indices, **kwargs)
    def select_trait(self, indices: ArrayLike, **kwargs: dict) -> TraitMatrix:
        return super().select_trait(indices, **kwargs)
    def sort(self, keys: tuple | numpy.ndarray, axis: int, **kwargs: dict) -> None:
        return super().sort(keys, axis, **kwargs)
    def sort_taxa(self, keys: tuple | numpy.ndarray, **kwargs: dict) -> None:
        return super().sort_taxa(keys, **kwargs)
    def sort_trait(self, keys: tuple | numpy.ndarray, **kwargs: dict) -> None:
        return super().sort_trait(keys, **kwargs)
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
    def trait(self) -> object:
        return super().trait
    @trait.setter
    def trait(self, value: object) -> None:
        super().trait = value
    @property
    def trait_axis(self) -> object:
        return super().trait_axis
    @trait_axis.setter
    def trait_axis(self, value: object) -> None:
        super().trait_axis = value
    def unscale(self) -> numpy.ndarray:
        return super().unscale()
    @property
    def location(self) -> object:
        return super().location
    @location.setter
    def location(self, value: object) -> None:
        super().location = value
    @property
    def scale(self) -> object:
        return super().scale
    @scale.setter
    def scale(self, value: object) -> None:
        super().scale = value
    def targmax(self) -> numpy.ndarray:
        return super().targmax()
    def targmin(self) -> numpy.ndarray:
        return super().targmin()
    def tmax(self, unscale: bool) -> numpy.ndarray:
        return super().tmax(unscale)
    def tmean(self, unscale: bool) -> numpy.ndarray:
        return super().tmean(unscale)
    def tmin(self, unscale: bool) -> numpy.ndarray:
        return super().tmin(unscale)
    def tstd(self, unscale: bool) -> numpy.ndarray:
        return super().tstd(unscale)
    def tvar(self, unscale: bool) -> numpy.ndarray:
        return super().tvar(unscale)
    def trange(self, unscale: bool) -> numpy.ndarray:
        return super().trange(unscale)
    @property
    def trait(self) -> object:
        return super().trait
    @trait.setter
    def trait(self, value: object) -> None:
        super().trait = value
    def to_csv(self, filename: str, **kwargs: dict) -> None:
        return super().to_csv(filename, **kwargs)
    def to_pandas(self, **kwargs: dict) -> pandas.DataFrame:
        return super().to_pandas(**kwargs)
    def ungroup(self, axis: int, **kwargs: dict) -> None:
        return super().ungroup(axis, **kwargs)
    def ungroup_taxa(self, **kwargs: dict) -> None:
        return super().ungroup_taxa(**kwargs)
    @classmethod
    def from_csv(cls, filename: str, **kwargs: dict) -> BreedingValueMatrix:
        return super().from_csv(filename, **kwargs)
    @classmethod
    def from_pandas(cls, df: pandas.DataFrame, **kwargs: dict) -> BreedingValueMatrix:
        return super().from_pandas(df, **kwargs)
    @classmethod
    def from_gmod(cls, gmod: GenomicModel, pgmat: PhasedGenotypeMatrix, nprogeny: Integral | numpy.ndarray, nrep: Integral | numpy.ndarray, **kwargs: dict) -> ExpectedMaximumBreedingValueMatrix:
        return super().from_gmod(gmod, pgmat, nprogeny, nrep, **kwargs)

