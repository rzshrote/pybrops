from typing import Optional

import pandas
from pybrops.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix


class DummyBreedingValueProtocol(BreedingValueProtocol):
    def estimate(self, ptobj: pandas.DataFrame, gtobj: GenotypeMatrix, miscout: dict | None, **kwargs: dict) -> BreedingValueMatrix:
        return super().estimate(ptobj, gtobj, miscout, **kwargs)
