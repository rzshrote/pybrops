"""
Module containing general error checking routines.
"""
# order dependent import

# generic
from .error_generic_python import *
from .error_generic_numpy import *

# python errors
from .error_attr_python import *
from .error_io_python import *
from .error_type_python import *
from .error_value_python import *

# h5py errors
from .error_io_h5py import *

# numpy errors
# from .error_attr_numpy import *
from .error_type_numpy import *
from .error_value_numpy import *

# pandas errors
from .error_type_pandas import *
