# 3rd party libraries

# our libraries
from . import Algorithm

class HillClimber(Algorithm):
    """docstring for HillClimber."""

    ############################################################################
    ######################### Reserved object methods ##########################
    ############################################################################
    def __init__(self, name = "Hill-Climber"):
        super(HillClimber, self).__init__(name)
