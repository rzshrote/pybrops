from numbers import Real
from typing import Union
from numpy import ndarray
import pandas
from pybrops.breed.prot.pt.PhenotypingProtocol import PhenotypingProtocol
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class DummyPhenotypingProtocol(PhenotypingProtocol):
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
    def phenotype(self, pgmat: PhasedGenotypeMatrix, miscout: dict, **kwargs: dict) -> pandas.DataFrame:
        return super().phenotype(pgmat, miscout, **kwargs)
    def set_H2(self, H2: Real | ndarray, pgmat: PhasedGenotypeMatrix, **kwargs: dict) -> None:
        return super().set_H2(H2, pgmat, **kwargs)
    def set_h2(self, h2: Real | ndarray, pgmat: PhasedGenotypeMatrix, **kwargs: dict) -> None:
        return super().set_h2(h2, pgmat, **kwargs)
    