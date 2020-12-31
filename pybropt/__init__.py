# order dependent import of submodules

# utils always first
from . import util

# core before all other submodules
from . import core

# order dependent submodules
from . import popgen
from . import model
from . import breed
