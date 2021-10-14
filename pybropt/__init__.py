# order dependent import of submodules

# core before all other submodules
from . import core

# algorithms
from . import algo

# order dependent submodules
from . import popgen
from . import model
from . import breed
