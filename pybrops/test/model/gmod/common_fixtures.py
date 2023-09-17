from typing import Optional, Union
from numpy import ndarray, dtype
import pandas
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.model.gmod.AdditiveDominanceEpistaticLinearGenomicModel import AdditiveDominanceEpistaticLinearGenomicModel
from pybrops.model.gmod.AdditiveDominanceLinearGenomicModel import AdditiveDominanceLinearGenomicModel
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.CoancestryLinearGenomicModel import CoancestryLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.gmod.NonlinearGenomicModel import NonlinearGenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.model.gmod.LinearGenomicModel import LinearGenomicModel

class DummyGenomicModel(GenomicModel):
    def __init__(self) -> None:
        pass
    def __copy__(self) -> GenomicModel:
        return super().__copy__()
    def __deepcopy__(self, memo: dict) -> GenomicModel:
        return super().__deepcopy__(memo)
    @property
    def model_name(self) -> object:
        return super().model_name
    @model_name.setter
    def model_name(self, value: object) -> None:
        super().model_name = value
    @property
    def params(self) -> object:
        return super().params
    @params.setter
    def params(self, value: object) -> None:
        super().params = value
    @property
    def trait(self) -> object:
        return super().trait
    @trait.setter
    def trait(self, value: object) -> None:
        super().trait = value
    @property
    def ntrait(self) -> object:
        return super().ntrait
    @ntrait.setter
    def ntrait(self, value: object) -> None:
        super().ntrait = value
    def fit_numpy(self, Y: ndarray, X: ndarray, Z: ndarray, **kwargs: dict) -> None:
        return super().fit_numpy(Y, X, Z, **kwargs)
    def fit(self, ptobj: BreedingValueMatrix | pandas.DataFrame | ndarray, cvobj: ndarray, gtobj: GenotypeMatrix | ndarray, **kwargs: dict) -> None:
        return super().fit(ptobj, cvobj, gtobj, **kwargs)
    def predict_numpy(self, X: ndarray, Z: ndarray, **kwargs: dict) -> ndarray:
        return super().predict_numpy(X, Z, **kwargs)
    def predict(self, cvobj: ndarray, gtobj: GenotypeMatrix, **kwargs: dict) -> BreedingValueMatrix:
        return super().predict(cvobj, gtobj, **kwargs)
    def score_numpy(self, Y: ndarray, X: ndarray, Z: ndarray, **kwargs: dict) -> ndarray:
        return super().score_numpy(Y, X, Z, **kwargs)
    def score(self, ptobj: BreedingValueMatrix | pandas.DataFrame, cvobj: object, gtobj: GenotypeMatrix, **kwargs: dict) -> ndarray:
        return super().score(ptobj, cvobj, gtobj, **kwargs)
    def gebv_numpy(self, Z: ndarray, **kwargs: dict) -> ndarray:
        return super().gebv_numpy(Z, **kwargs)
    def gebv(self, gtobj: GenotypeMatrix, **kwargs: dict) -> BreedingValueMatrix:
        return super().gebv(gtobj, **kwargs)
    def var_G_numpy(self, Z: ndarray, **kwargs: dict) -> ndarray:
        return super().var_G_numpy(Z, **kwargs)
    def var_G(self, gtobj: GenotypeMatrix, **kwargs: dict) -> ndarray:
        return super().var_G(gtobj, **kwargs)
    def var_A_numpy(self, Z: ndarray, **kwargs: dict) -> ndarray:
        return super().var_A_numpy(Z, **kwargs)
    def var_A(self, gtobj: GenotypeMatrix, **kwargs: dict) -> ndarray:
        return super().var_A(gtobj, **kwargs)
    def var_a_numpy(self, p: ndarray, ploidy: int, **kwargs: dict) -> ndarray:
        return super().var_a_numpy(p, ploidy, **kwargs)
    def var_a(self, gtobj: GenotypeMatrix | ndarray, ploidy: int, **kwargs: dict) -> ndarray:
        return super().var_a(gtobj, ploidy, **kwargs)
    def bulmer_numpy(self, Z: ndarray, p: ndarray, ploidy: int, **kwargs: dict) -> ndarray:
        return super().bulmer_numpy(Z, p, ploidy, **kwargs)
    def bulmer(self, gtobj: GenotypeMatrix, ploidy: int, **kwargs: dict) -> ndarray:
        return super().bulmer(gtobj, ploidy, **kwargs)
    def usl_numpy(self, p: ndarray, ploidy: int, descale: bool, **kwargs: dict) -> ndarray:
        return super().usl_numpy(p, ploidy, descale, **kwargs)
    def usl(self, gtobj: GenotypeMatrix, ploidy: int, descale: bool, **kwargs: dict) -> ndarray:
        return super().usl(gtobj, ploidy, descale, **kwargs)
    def lsl_numpy(self, p: ndarray, ploidy: int, descale: bool, **kwargs: dict) -> ndarray:
        return super().lsl_numpy(p, ploidy, descale, **kwargs)
    def lsl(self, gtobj: GenotypeMatrix, ploidy: int, descale: bool, **kwargs: dict) -> ndarray:
        return super().lsl(gtobj, ploidy, descale, **kwargs)
    def facount(self, gmat: GenotypeMatrix, dtype: dtype | None, **kwargs: dict) -> ndarray:
        return super().facount(gmat, dtype, **kwargs)
    def fafreq(self, gmat: GenotypeMatrix, dtype: ndarray | None, **kwargs: dict) -> ndarray:
        return super().fafreq(gmat, dtype, **kwargs)
    def faavail(self, gmat: GenotypeMatrix, dtype: dtype | None = None, **kwargs: dict) -> ndarray:
        return super().faavail(gmat, dtype, **kwargs)
    def fafixed(self, gmat: GenotypeMatrix, dtype: ndarray | None, **kwargs: dict) -> ndarray:
        return super().fafixed(gmat, dtype, **kwargs)
    def dacount(self, gmat: GenotypeMatrix, dtype: ndarray | None, **kwargs: dict) -> ndarray:
        return super().dacount(gmat, dtype, **kwargs)
    def dafreq(self, gmat: GenotypeMatrix, dtype: ndarray | None, **kwargs: dict) -> ndarray:
        return super().dafreq(gmat, dtype, **kwargs)
    def daavail(self, gmat: GenotypeMatrix, dtype: dtype | None = None, **kwargs: dict) -> ndarray:
        return super().daavail(gmat, dtype, **kwargs)
    def dafixed(self, gmat: GenotypeMatrix, dtype: ndarray | None, **kwargs: dict) -> ndarray:
        return super().dafixed(gmat, dtype, **kwargs)
    @classmethod
    def from_hdf5(cls, filename: str, groupname: str | None) -> HDF5InputOutput:
        return super().from_hdf5(filename, groupname)
    def to_hdf5(self, filename: str, groupname: str | None) -> None:
        return super().to_hdf5(filename, groupname)

class DummyLinearGenomicModel(DummyGenomicModel,LinearGenomicModel):
    @property
    def beta(self) -> object:
        return super().beta
    @beta.setter
    def beta(self, value: object) -> None:
        super().beta = value
    @property
    def u(self) -> object:
        return super().u
    @u.setter
    def u(self, value: object) -> None:
        super().u = value

class DummyNonlinearGenomicModel(DummyGenomicModel,NonlinearGenomicModel):
    pass

class DummyAdditiveLinearGenomicModel(DummyLinearGenomicModel,AdditiveLinearGenomicModel):
    @property
    def u_misc(self) -> object:
        return super().u_misc
    @u_misc.setter
    def u_misc(self, value: object) -> None:
        super().u_misc = value
    @property
    def u_a(self) -> object:
        return super().u_a
    @u_a.setter
    def u_a(self, value: object) -> None:
        super().u_a = value

class DummyAdditiveDominanceLinearGenomicModel(DummyAdditiveLinearGenomicModel,AdditiveDominanceLinearGenomicModel):
    @property
    def u_d(self) -> object:
        return super().u_d
    @u_d.setter
    def u_d(self, value: object) -> None:
        super().u_d = value

class DummyAdditiveDominanceEpistaticLinearGenomicModel(DummyAdditiveDominanceLinearGenomicModel,AdditiveDominanceEpistaticLinearGenomicModel):
    @property
    def u_i(self) -> object:
        return super().u_i
    @u_i.setter
    def u_i(self, value: object) -> None:
        super().u_i = value

class DummyCoancestryLinearGenomicModel(DummyLinearGenomicModel,CoancestryLinearGenomicModel):
    @property
    def u_misc(self) -> object:
        return super().u_misc
    @u_misc.setter
    def u_misc(self, value: object) -> None:
        super().u_misc = value
    @property
    def u_c(self) -> object:
        return super().u_c
    @u_c.setter
    def u_c(self, value: object) -> None:
        super().u_c = value
