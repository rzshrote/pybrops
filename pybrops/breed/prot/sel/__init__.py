"""
Module containing selection protocols.
"""

__all__ = [
    "cfg",
    "prob",
    "sampling",
    "targetfn",
    "transfn",
    "weightfn",
    "SelectionProtocol",
    "UnconstrainedSelectionProtocol",
    "EstimatedBreedingValueSelection",
    "ExpectedMaximumBreedingValueSelection",
    "FamilyEstimatedBreedingValueSelection",
    "GeneralizedWeightedGenomicEstimatedBreedingValueSelection",
    "GenomicEstimatedBreedingValueSelection",
    "GenotypeBuilderSelection",
    "L1NormGenomicSelection",
    "L2NormGenomicSelection",
    "MeanExpectedHeterozygositySelection",
    "MeanGenomicRelationshipSelection",
    "MultiObjectiveGenomicSelection",
    "OptimalContributionSelection",
    "OptimalHaploidValueSelection",
    "OptimalPopulationValueSelection",
    "RandomSelection",
    "UsefulnessCriterionSelection",
    "WeightedGenomicSelection"
]

# order dependent imports

# import submodules
from pybrops.breed.prot.sel import cfg
from pybrops.breed.prot.sel import prob

# import utilities
from pybrops.breed.prot.sel import targetfn
from pybrops.breed.prot.sel import transfn
from pybrops.breed.prot.sel import weightfn

# abstract classes
from pybrops.breed.prot.sel import SelectionProtocol
from pybrops.breed.prot.sel import UnconstrainedSelectionProtocol

# concrete classes
from pybrops.breed.prot.sel import EstimatedBreedingValueSelection
from pybrops.breed.prot.sel import ExpectedMaximumBreedingValueSelection
from pybrops.breed.prot.sel import FamilyEstimatedBreedingValueSelection
from pybrops.breed.prot.sel import GeneralizedWeightedGenomicEstimatedBreedingValueSelection
from pybrops.breed.prot.sel import GenomicEstimatedBreedingValueSelection
from pybrops.breed.prot.sel import GenotypeBuilderSelection
from pybrops.breed.prot.sel import L1NormGenomicSelection
from pybrops.breed.prot.sel import L2NormGenomicSelection
from pybrops.breed.prot.sel import MeanExpectedHeterozygositySelection
from pybrops.breed.prot.sel import MeanGenomicRelationshipSelection
from pybrops.breed.prot.sel import MultiObjectiveGenomicSelection
from pybrops.breed.prot.sel import OptimalContributionSelection
from pybrops.breed.prot.sel import OptimalHaploidValueSelection
from pybrops.breed.prot.sel import OptimalPopulationValueSelection
from pybrops.breed.prot.sel import RandomSelection
from pybrops.breed.prot.sel import UsefulnessCriterionSelection
from pybrops.breed.prot.sel import WeightedGenomicSelection
