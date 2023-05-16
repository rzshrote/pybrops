"""
Module containing general error checking routines.
"""

__all__ = [
    "error_generic_python",
    "error_generic_numpy",
    "error_attr_python",
    "error_io_python",
    "error_type_python",
    "error_value_python",
    "error_io_h5py",
    "error_type_numpy",
    "error_value_numpy",
    "error_type_pandas"
]

# order dependent import

# generic
from pybrops.core.error import error_generic_python
from pybrops.core.error import error_generic_numpy

# python errors
from pybrops.core.error import error_attr_python
from pybrops.core.error import error_io_python
from pybrops.core.error import error_type_python
from pybrops.core.error import error_value_python

# h5py errors
from pybrops.core.error import error_io_h5py

# numpy errors
from pybrops.core.error import error_type_numpy
from pybrops.core.error import error_value_numpy

# pandas errors
from pybrops.core.error import error_type_pandas
