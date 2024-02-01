from typing import Optional, Union
from numpy import ndarray
import pandas
from pybrops.popgen.gmap.GeneticMap import GeneticMap
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction


class DummyGeneticMapFunction(GeneticMapFunction):
    def __init__(self) -> None:
        pass
    def invmapfn(self, r: ndarray) -> ndarray:
        return super().invmapfn(r)
    def mapfn(self, d: ndarray) -> ndarray:
        return super().mapfn(d)
    def rprob1g(self, gmap: GeneticMap, vrnt_chrgrp: ndarray, vrnt_genpos: ndarray) -> ndarray:
        return super().rprob1g(gmap, vrnt_chrgrp, vrnt_genpos)
    def rprob1p(self, gmap: GeneticMap, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray) -> ndarray:
        return super().rprob1p(gmap, vrnt_chrgrp, vrnt_phypos)
    def rprob2g(self, gmap: GeneticMap, vrnt_chrgrp: ndarray, vrnt_genpos: ndarray) -> ndarray:
        return super().rprob2g(gmap, vrnt_chrgrp, vrnt_genpos)
    def rprob2p(self, gmap: GeneticMap, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray) -> ndarray:
        return super().rprob2p(gmap, vrnt_chrgrp, vrnt_phypos)

class DummyGeneticMap(GeneticMap):
    def __init__(self) -> None:
        pass
    def __len__(self):
        return super().__len__()
    def build_spline(self, kind: str, fill_value: object, **kwargs: dict) -> None:
        return super().build_spline(kind, fill_value, **kwargs)
    def congruence(self) -> ndarray:
        return super().congruence()
    def gdist1g(self, vrnt_chrgrp: ndarray, vrnt_genpos: ndarray, ast: int | None, asp: int | None) -> ndarray:
        return super().gdist1g(vrnt_chrgrp, vrnt_genpos, ast, asp)
    def gdist1p(self, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray, ast: int | None, asp: int | None) -> ndarray:
        return super().gdist1p(vrnt_chrgrp, vrnt_phypos, ast, asp)
    def gdist2g(self, vrnt_chrgrp: ndarray, vrnt_genpos: ndarray, rst: int | None, rsp: int | None, cst: int | None, csp: int | None) -> ndarray:
        return super().gdist2g(vrnt_chrgrp, vrnt_genpos, rst, rsp, cst, csp)
    def gdist2p(self, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray, rst: int | None, rsp: int | None, cst: int | None, csp: int | None) -> ndarray:
        return super().gdist2p(vrnt_chrgrp, vrnt_phypos, rst, rsp, cst, csp)
    def group(self, **kwargs: dict) -> None:
        return super().group(**kwargs)
    def has_spline(self) -> bool:
        return super().has_spline()
    def interp_genpos(self, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray) -> ndarray:
        return super().interp_genpos(vrnt_chrgrp, vrnt_phypos)
    def interp_gmap(self, vrnt_chrgrp: ndarray, vrnt_phypos: ndarray, **kwargs: dict):
        return super().interp_gmap(vrnt_chrgrp, vrnt_phypos, **kwargs)
    def is_congruent(self) -> bool:
        return super().is_congruent()
    def is_grouped(self, **kwargs: dict) -> bool:
        return super().is_grouped(**kwargs)
    def lexsort(self, keys, **kwargs: dict) -> ndarray:
        return super().lexsort(keys, **kwargs)
    def prune(self, nt: int, M: float, **kwargs: dict) -> None:
        return super().prune(nt, M, **kwargs)
    def remove(self, indices, **kwargs: dict) -> None:
        return super().remove(indices, **kwargs)
    def remove_discrepancies(self) -> None:
        return super().remove_discrepancies()
    def reorder(self, indices, **kwargs: dict) -> None:
        return super().reorder(indices, **kwargs)
    def select(self, indices, **kwargs: dict) -> None:
        return super().select(indices, **kwargs)
    def sort(self, keys, **kwargs: dict) -> None:
        return super().sort(keys, **kwargs)
    @property
    def spline(self) -> object:
        return super().spline
    @spline.setter
    def spline(self, value: object) -> None:
        super().spline = value
    @property
    def spline_kind(self) -> object:
        return super().spline_kind
    @spline_kind.setter
    def spline_kind(self, value: object) -> None:
        super().spline_kind = value
    @property
    def spline_fill_value(self) -> object:
        return super().spline_fill_value
    @spline_fill_value.setter
    def spline_fill_value(self, value: object) -> None:
        super().spline_fill_value = value
    @property
    def vrnt_chrgrp(self) -> object:
        return super().vrnt_chrgrp
    @vrnt_chrgrp.setter
    def vrnt_chrgrp(self, value: object) -> None:
        super().vrnt_chrgrp = value
    @property
    def vrnt_phypos(self) -> object:
        return super().vrnt_phypos
    @vrnt_phypos.setter
    def vrnt_phypos(self, value: object) -> None:
        super().vrnt_phypos = value
    @property
    def vrnt_genpos(self) -> object:
        return super().vrnt_genpos
    @vrnt_genpos.setter
    def vrnt_genpos(self, value: object) -> None:
        super().vrnt_genpos = value
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
    def vrnt_chrgrp_len(self) -> object:
        return super().vrnt_chrgrp_len
    @vrnt_chrgrp_len.setter
    def vrnt_chrgrp_len(self, value: object) -> None:
        super().vrnt_chrgrp_len = value
    def select(self, indices, **kwargs: dict) -> None:
        return super().select(indices, **kwargs)
    def sort(self, keys, **kwargs: dict) -> None:
        return super().sort(keys, **kwargs)
    def to_csv(self, filename: str, sep: str, header: bool, index: bool | int, **kwargs: dict) -> None:
        return super().to_csv(filename, sep, header, index, **kwargs)
    def to_pandas(self) -> pandas.DataFrame:
        return super().to_pandas()
    def __copy__(self) -> GeneticMap:
        return super().__copy__()
    def __deepcopy__(self, memo: dict) -> GeneticMap:
        return super().__deepcopy__(memo)
    def copy(self) -> GeneticMap:
        return super().copy()
    def deepcopy(self, memo: dict) -> GeneticMap:
        return super().deepcopy(memo)
    @classmethod
    def from_csv(cls, filename: str, **kwargs: dict) -> GeneticMap:
        return super().from_csv(filename, **kwargs)
    @classmethod
    def from_pandas(cls, df: pandas.DataFrame, **kwargs: dict) -> GeneticMap:
        return super().from_pandas(df, **kwargs)
    @property
    def nvrnt(self) -> object:
        return super().nvrnt
    def ungroup(self, **kwargs: dict) -> None:
        return super().ungroup(**kwargs)
