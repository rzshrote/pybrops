from numbers import Real
from typing import Sequence
from numpy import ndarray
from numpy.typing import ArrayLike, DTypeLike
import pandas
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.core.mat import TaxaMatrix
from pybrops.core.mat.Matrix import Matrix
from pybrops.popgen.cmat.CoancestryMatrix import CoancestryMatrix
from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix


class DummyCoancestryMatrix(CoancestryMatrix):
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
    def adjoin(self, values: Matrix | ndarray, axis: int, **kwargs: dict) -> Matrix:
        return super().adjoin(values, axis, **kwargs)
    def adjoin_taxa(self, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> TaxaMatrix:
        return super().adjoin_taxa(values, taxa, taxa_grp, **kwargs)
    def append(self, values: Matrix | ndarray, axis: int, **kwargs: dict) -> None:
        return super().append(values, axis, **kwargs)
    def append_taxa(self, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> None:
        return super().append_taxa(values, taxa, taxa_grp, **kwargs)
    @classmethod
    def concat(cls, mats: Sequence, axis: int, **kwargs: dict) -> Matrix:
        return super().concat(mats, axis, **kwargs)
    @classmethod
    def concat_taxa(cls, mats: Sequence, **kwargs: dict) -> TaxaMatrix:
        return super().concat_taxa(mats, **kwargs)
    def copy(self) -> Matrix:
        return super().copy()
    def deepcopy(self, memo: dict) -> Matrix:
        return super().deepcopy(memo)
    def delete(self, obj: int | slice | Sequence, axis: int, **kwargs: dict) -> Matrix:
        return super().delete(obj, axis, **kwargs)
    def delete_taxa(self, obj: int | slice | Sequence, **kwargs: dict) -> TaxaMatrix:
        return super().delete_taxa(obj, **kwargs)
    @classmethod
    def from_hdf5(cls, filename: str, groupname: str | None) -> HDF5InputOutput:
        return super().from_hdf5(filename, groupname)
    def group(self, axis: int, **kwargs: dict) -> None:
        return super().group(axis, **kwargs)
    def group_taxa(self, **kwargs: dict) -> None:
        return super().group_taxa(**kwargs)
    def incorp(self, obj: int | slice | Sequence, values: Matrix | ndarray, axis: int, **kwargs) -> None:
        return super().incorp(obj, values, axis, **kwargs)
    def incorp_taxa(self, obj: int | slice | Sequence, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> None:
        return super().incorp_taxa(obj, values, taxa, taxa_grp, **kwargs)
    def insert(self, obj: int | slice | Sequence, values: ArrayLike, axis: int, **kwargs: dict) -> Matrix:
        return super().insert(obj, values, axis, **kwargs)
    def insert_taxa(self, obj: int | slice | Sequence, values: Matrix | ndarray, taxa: ndarray, taxa_grp: ndarray, **kwargs: dict) -> TaxaMatrix:
        return super().insert_taxa(obj, values, taxa, taxa_grp, **kwargs)
    def is_grouped(self, axis: int, **kwargs: dict) -> bool:
        return super().is_grouped(axis, **kwargs)
    def is_grouped_taxa(self, **kwargs: dict) -> bool:
        return super().is_grouped_taxa(**kwargs)
    def is_square(self) -> bool:
        return super().is_square()
    def lexsort(self, keys: tuple | ndarray, axis: int, **kwargs: dict) -> ndarray:
        return super().lexsort(keys, axis, **kwargs)
    def lexsort_taxa(self, keys: tuple | ndarray, **kwargs: dict) -> ndarray:
        return super().lexsort_taxa(keys, **kwargs)
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
    def nsquare(self) -> object:
        return super().nsquare
    @property
    def ntaxa(self) -> object:
        return super().ntaxa
    @ntaxa.setter
    def ntaxa(self, value: object) -> None:
        super().ntaxa = value
    def remove(self, obj: int | slice | Sequence, axis: int, **kwargs: dict) -> None:
        return super().remove(obj, axis, **kwargs)
    def remove_taxa(self, obj: int | slice | Sequence, **kwargs: dict) -> None:
        return super().remove_taxa(obj, **kwargs)
    def reorder(self, indices: ndarray | Sequence, axis: int, **kwargs: dict) -> None:
        return super().reorder(indices, axis, **kwargs)
    def reorder_taxa(self, indices: ndarray | Sequence, **kwargs: dict) -> None:
        return super().reorder_taxa(indices, **kwargs)
    def select(self, indices: ArrayLike, axis: int, **kwargs: dict) -> Matrix:
        return super().select(indices, axis, **kwargs)
    def select_taxa(self, indices: ArrayLike, **kwargs: dict) -> TaxaMatrix:
        return super().select_taxa(indices, **kwargs)
    def sort(self, keys: tuple | ndarray, axis: int, **kwargs: dict) -> None:
        return super().sort(keys, axis, **kwargs)
    def sort_taxa(self, keys: tuple | ndarray, **kwargs: dict) -> None:
        return super().sort_taxa(keys, **kwargs)
    @property
    def square_axes(self) -> object:
        return super().square_axes
    @square_axes.setter
    def square_axes(self, value: object) -> None:
        super().square_axes = value
    @property
    def square_axes_len(self) -> object:
        return super().square_axes_len
    @square_axes_len.setter
    def square_axes_len(self, value: object) -> None:
        super().square_axes_len = value
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
    def apply_jitter(self, eigvaltol: float, minjitter: float, maxjitter: float, nattempt: int) -> bool:
        return super().apply_jitter(eigvaltol, minjitter, maxjitter, nattempt)
    def coancestry(self, *args: tuple, **kwargs: dict) -> Real:
        return super().coancestry(*args, **kwargs)
    @classmethod
    def from_gmat(cls, gmat: GenotypeMatrix, **kwargs: dict) -> CoancestryMatrix:
        return super().from_gmat(gmat, **kwargs)
    def inverse(self, format: str) -> ndarray:
        return super().inverse(format)
    def is_positive_semidefinite(self, eigvaltol: float) -> bool:
        return super().is_positive_semidefinite(eigvaltol)
    def kinship(self, *args: tuple, **kwargs: dict) -> Real:
        return super().kinship(*args, **kwargs)
    def mat_asformat(self, format: str) -> ndarray:
        return super().mat_asformat(format)
    def max(self, format: str, axis: int | tuple | None) -> Real | ndarray:
        return super().max(format, axis)
    def max_inbreeding(self, format: str) -> Real:
        return super().max_inbreeding(format)
    def mean(self, format: str, axis: int | tuple | None, dtype: DTypeLike) -> Real:
        return super().mean(format, axis, dtype)
    def min(self, format: str, axis: int | tuple | None) -> Real | ndarray:
        return super().min(format, axis)
    def min_inbreeding(self, format: str) -> Real:
        return super().min_inbreeding(format)
    def to_pandas(self, **kwargs: dict) -> pandas.DataFrame:
        return super().to_pandas(**kwargs)
    def to_csv(self, filename: str, **kwargs: dict) -> None:
        return super().to_csv(filename, **kwargs)
    @classmethod
    def from_pandas(cls, df: pandas.DataFrame, **kwargs: dict) -> CoancestryMatrix:
        return super().from_pandas(df, **kwargs)
    @classmethod
    def from_csv(cls, filename: str, **kwargs: dict) -> CoancestryMatrix:
        return super().from_csv(filename, **kwargs)
    def ungroup(self, axis: int, **kwargs: dict) -> None:
        return super().ungroup(axis, **kwargs)
    def ungroup_taxa(self, **kwargs: dict) -> None:
        return super().ungroup_taxa(**kwargs)

class DummyDenseCoancestryMatrix(DenseCoancestryMatrix):
    @classmethod
    def from_gmat(cls, gmat: GenotypeMatrix, **kwargs: dict) -> CoancestryMatrix:
        return super().from_gmat(gmat, **kwargs)
    