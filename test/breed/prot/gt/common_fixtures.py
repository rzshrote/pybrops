from typing import Optional
from pybrops.breed.prot.gt.GenotypingProtocol import GenotypingProtocol
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class DummyGenotypingProtocol(GenotypingProtocol):
    def genotype(self, pgmat: PhasedGenotypeMatrix, miscout: dict | None, **kwargs: dict) -> GenotypeMatrix:
        return super().genotype(pgmat, miscout, **kwargs)
