"""
Module containing mating simulation protocols.
"""

__all__ = [
    "util",
    "MatingProtocol",
    "SelfCross",
    "TwoWayCross",
    "TwoWayDHCross",
    "ThreeWayCross",
    "ThreeWayDHCross",
    "FourWayCross",
    "FourWayDHCross",
]

# import libraries; order dependent

# import helper functions
from pybrops.breed.prot.mate import util

# abstract classes
from pybrops.breed.prot.mate import MatingProtocol

# concrete classes
from pybrops.breed.prot.mate import SelfCross
from pybrops.breed.prot.mate import TwoWayCross
from pybrops.breed.prot.mate import TwoWayDHCross
from pybrops.breed.prot.mate import ThreeWayCross
from pybrops.breed.prot.mate import ThreeWayDHCross
from pybrops.breed.prot.mate import FourWayCross
from pybrops.breed.prot.mate import FourWayDHCross
