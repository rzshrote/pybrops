from numbers import Real
from pathlib import Path
from typing import Union
from h5py._hl.files import File
from numpy import ndarray
import pandas
from pybrops.breed.prot.pt.PhenotypingProtocol import PhenotypingProtocol
from pybrops.core.io.Copyable import Copyable
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class DummyPhenotypingProtocol(PhenotypingProtocol):
    def __copy__(self) -> Copyable:
        return super().__copy__()
    def __deepcopy__(self, memo: dict | None) -> Copyable:
        return super().__deepcopy__(memo)
    @property
    def gpmod(self) -> object:
        """gpmod."""
        return self._gpmod
    @gpmod.setter
    def gpmod(self, value: object) -> None:
        """Set gpmod."""
        self._gpmod = value
    @property
    def var_err(self) -> object:
        """var_err."""
        return self._var_err
    @var_err.setter
    def var_err(self, value: object) -> None:
        """Set var_err."""
        self._var_err = value
    def copy(self) -> Copyable:
        return super().copy()
    def deepcopy(self, memo: dict | None) -> Copyable:
        return super().deepcopy(memo)
    def phenotype(self, pgmat: PhasedGenotypeMatrix, miscout: dict, **kwargs: dict) -> pandas.DataFrame:
        return super().phenotype(pgmat, miscout, **kwargs)
    def set_H2(self, H2: Real | ndarray, pgmat: PhasedGenotypeMatrix, **kwargs: dict) -> None:
        return super().set_H2(H2, pgmat, **kwargs)
    def set_h2(self, h2: Real | ndarray, pgmat: PhasedGenotypeMatrix, **kwargs: dict) -> None:
        return super().set_h2(h2, pgmat, **kwargs)
    def to_hdf5(self, filename: str | Path | File, groupname: str | None, overwrite: bool) -> None:
        return super().to_hdf5(filename, groupname, overwrite)
    @classmethod
    def from_hdf5(cls, filename: str | Path | File, groupname: str | None) -> HDF5InputOutput:
        return super().from_hdf5(filename, groupname)
    