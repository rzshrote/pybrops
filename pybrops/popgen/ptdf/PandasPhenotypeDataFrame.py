"""
Module implementing a phenotype dataframe using Pandas DataFrames and its
associated error checking routines.
"""

from typing import Any
import numpy
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_python import check_len
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame
from pybrops.core.df.PandasDataFrame import PandasDataFrame
from pybrops.core.error.error_type_pandas import check_is_pandas_df

class PandasPhenotypeDataFrame(PandasDataFrame,PhenotypeDataFrame):
    """
    A concrete class for data frame objects utilizing Pandas DataFrames as a
    storage container.
    """

    ########################## Special Object Methods ##########################
    def __init__(self, data, col_atype, col_aefct, **kwargs: dict):
        super(PandasPhenotypeDataFrame, self).__init__(
            data = data,
            **kwargs
        )
        self.col_atype = col_atype
        self.col_aefct = col_aefct

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################## Column attributes ###################
    @property
    def col_ctype(self) -> numpy.ndarray:
        """Column types used for classifying variables."""
        return self._col_ctype
    @col_ctype.setter
    def col_ctype(self, value: numpy.ndarray) -> None:
        """Set column types"""
        check_is_ndarray(value, "col_ctype")
        check_len(value, "col_ctype", self.ncol)
        options = ["response","predictor",None]
        if any(e not in options for e in value):
            raise ValueError("elements in 'col_ctype' must be in {0}".format(options))
        self._col_ctype = value
    @col_ctype.deleter
    def col_ctype(self) -> None:
        """Delete column types"""
        del self._col_ctype

    @property
    def col_atype(self) -> numpy.ndarray:
        """Analysis variable type array."""
        return self._col_atype
    @col_atype.setter
    def col_atype(self, value: numpy.ndarray) -> None:
        """Set analysis variable type array"""
        check_is_ndarray(value, "col_atype")
        check_len(value, "col_atype", self.ncol)
        options = ['bool', 'complex', 'double', 'int', 'raw', 'str',
            'factor(bool)', 'factor(complex)', 'factor(double)',
            'factor(int)', 'factor(str)', None
        ]
        if any(e not in options for e in value):
            raise ValueError("elements in 'col_atype' must be in {0}".format(options))
        self._col_atype = value
    @col_atype.deleter
    def col_atype(self) -> None:
        """Delete analysis variable type array"""
        del self._col_atype

    @property
    def col_aefct(self) -> numpy.ndarray:
        """Analysis variable effect type {'response', 'fixed', 'random', None} array."""
        return self._col_aefct
    @col_aefct.setter
    def col_aefct(self, value: numpy.ndarray) -> None:
        """Set analysis variable effect type array"""
        check_is_ndarray(value, "col_aefct")
        check_len(value, "col_aefct", self.ncol)
        col_aefct_options = ['response', 'fixed', 'random', None]
        if any(e not in col_aefct_options for e in value):
            raise ValueError("elements in 'col_aefct' must be in {0}".format(col_aefct_options))
        self._col_aefct = value
    @col_aefct.deleter
    def col_aefct(self) -> None:
        """Delete analysis variable effect type array"""
        del self._col_aefct

    ############################## Object Methods ##############################
    def col_data(
            self, 
            index = None, 
            name = None, 
            ctype = None, 
            dtype = None, 
            atype = None, 
            aefct = None, 
            return_index = False, 
            return_name = False, 
            return_ctype = False, 
            return_dtype = False, 
            return_atype = False, 
            return_aefct = False, 
            **kwargs: dict
        ):
        """
        Get a column's (or columns') data from the dataframe.

        Parameters
        ----------
        index : int, None
            Integer index of the column to get.
        name : str, None
            Name of the column to get.
        ctype : str, None
            Classification type of the column to get.
        dtype : str, numpy.dtype, None
            Data type of the column to get.
        atype : str, None
            Analysis type of the column to get.
        aefct : str, None
            Analysis effect type of the column to get.
        return_index : boolean, default = False
            Whether to return the column index along with the column data.
        return_name : boolean, default = False
            Whether to return the column name along with the column data.
        return_ctype : boolean, default = False
            Whether to return the column type along with the column data.
        return_dtype : boolean, default = False
            Whether to return the column dtype along with the column data.
        return_atype : boolean, default = False
            Whether to return the column analysis type along with the column
            data.
        return_aefct : boolean, default = False
            Whether to return the column analysis effect type along with the
            column data.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray, tuple, None
            If no other returns are specified, return a numpy.ndarray for the
            column. Otherwise return a tuple of values corresponding to (in
            order of the function API) the requested data.
        """
        # construct a selection mask for column names
        mask = True
        if index is not None:
            tmp = numpy.empty(self.ncol, dtype = 'bool')
            tmp[index] = True
            mask = mask & tmp
        if name is not None:
            mask = mask & (self.col_name == name)
        if ctype is not None:
            mask = mask & (self.col_ctype == ctype)
        if dtype is not None:
            mask = mask & (self.col_dtype == dtype)
        if atype is not None:
            mask = mask & (self.col_atype == atype)
        if aefct is not None:
            mask = mask & (self.col_aefct == aefct)
        if mask is True:
            mask = False

        # get index of column equal to the provided name
        cix = numpy.flatnonzero(self.col_name == name)

        # process errors if key is not found or multiple keys are found
        if len(cix) == 0:
            err =  "PandasPhenotypeDataFrame does not contain a column with matches for: \n"
            err += "    name = {0}".format(name)
            err += "    ctype = {0}".format(ctype)
            err += "    dtype = {0}".format(dtype)
            raise KeyError(err)

        # extract column arrays as a numpy.ndarray
        out_arr = [self.df.iloc[:,ix].values for ix in cix]

        # construct extra output list
        out_extra = []

        # add values to extra output
        if return_index:
            out_extra.append(cix)
        if return_name:
            out_extra.append(self.col_name[cix])
        if return_ctype:
            out_extra.append(self.col_ctype[cix])
        if return_dtype:
            out_extra.append(self.col_dtype[cix])
        if return_atype:
            out_extra.append(self.col_atype[cix])
        if return_aefct:
            out_extra.append(self.col_aefct[cix])

        # construct output object (numpy.ndarray or tuple)
        out = out_arr if len(out_extra) == 0 else (out_arr, *out_extra)

        return out



################################################################################
################################## Utilities ###################################
################################################################################
def is_PandasPhenotypeDataFrame(v: object) -> bool:
    """
    Determine whether an object is a PandasPhenotypeDataFrame.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a PandasPhenotypeDataFrame object instance.
    """
    return isinstance(v, PandasPhenotypeDataFrame)

def check_is_PandasPhenotypeDataFrame(v: object, vname: str) -> None:
    """
    Check if object is of type PandasPhenotypeDataFrame. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, PandasPhenotypeDataFrame):
        raise TypeError("variable '{0}' must be a PandasPhenotypeDataFrame".format(vname))
