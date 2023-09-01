from typing import Optional
from pybrops.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame


class DummyBreedingValueProtocol(BreedingValueProtocol):
    def estimate(self, ptobj: PhenotypeDataFrame, gtobj: GenotypeMatrix, miscout: dict | None, **kwargs: dict) -> BreedingValueMatrix:
        return super().estimate(ptobj, gtobj, miscout, **kwargs)
