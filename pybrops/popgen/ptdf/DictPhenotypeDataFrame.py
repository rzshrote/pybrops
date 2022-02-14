import numpy
import copy

from pybrops.core.df.DictDataFrame import DictDataFrame
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

from pybrops.core.error import check_keys_in_dict
from pybrops.core.error import check_values_in_dict_all_type
from pybrops.core.error import check_len

class DictPhenotypeDataFrame(DictDataFrame,PhenotypeDataFrame):
    """Concrete class for phenotype dataframe objects."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, data, col_grp = None, col_analysis_type = None, col_analysis_effect = None, row_name = None, **kwargs):
        """
        Constructor for the concrete class DictPhenotypeDataFrame.

        Parameters
        ----------
        data : dict
        col_grp : dict
        col_analysis_type : dict
        col_analysis_effect : dict
        row_name : dict
        kwargs : dict
            Additional keyword arguments.
        """
        super(DictPhenotypeDataFrame, self).__init__(
            data = data,
            col_grp = col_grp,
            row_name = row_name,
            **kwargs
        )

        ### error checks and assignments (order dependent)
        self.col_analysis_type = col_analysis_type
        self.col_analysis_effect = col_analysis_effect

    def __copy__(self):
        """
        Make a shallow copy of the the dataframe.

        Returns
        -------
        out : DataFrame
        """
        return self.__class__(
            data = copy.copy(self.data),
            col_grp = copy.copy(self.col_grp),
            col_analysis_type = copy.copy(self.col_analysis_type),
            col_analysis_effect = copy.copy(self.col_analysis_effect),
            row_name = copy.copy(self.row_name)
        )

    def __deepcopy__(self, memo):
        """
        Make a deep copy of the dataframe.

        Parameters
        ----------
        memo : dict

        Returns
        -------
        out : DataFrame
        """
        return self.__class__(
            data = copy.deepcopy(self.data),
            col_grp = copy.deepcopy(self.col_grp),
            col_analysis_type = copy.deepcopy(self.col_analysis_type),
            col_analysis_effect = copy.deepcopy(self.col_analysis_effect),
            row_name = copy.deepcopy(self.row_name)
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def col_analysis_type():
        doc = "Analysis variable type array."
        def fget(self):
            """Get analysis types as numpy.ndarray or None if empty."""
            if hasattr(self, "_col_analysis_type") and isinstance(self._col_analysis_type, dict):
                return numpy.object_(list(self._col_analysis_type.values()))
            return None
        def fset(self, value):
            """Set analysis variable type array"""
            if isinstance(value, (list,tuple,numpy.ndarray)):   # if array_like
                check_len(value, "col_analysis_type", self.ncol)# check input length
                names = self.col_name                           # get column names
                grps = value                                    # get column grps
                value = dict(zip(names,grps))                   # construct dict of grps
            if value is None:                                   # if is None
                self._col_analysis_type = value                 # set to None
            elif isinstance(value, dict):                       # if is dict
                check_keys_in_dict(value, "col_analysis_type", *self._data.keys())
                check_values_in_dict_all_type(value, "col_analysis_type", (str,type(None)))
                options = ['bool', 'complex', 'double', 'int', 'raw', 'str',
                    'factor(bool)', 'factor(complex)', 'factor(double)',
                    'factor(int)', 'factor(str)', None
                ]
                if any(e not in options for e in value.values()):
                    raise ValueError("unsupported value: supported values are {0}".format(options))
                self._col_analysis_type = {k: value[k] for k in self._data.keys()}
            else:
                raise TypeError("unsupported type: supported types are dict, list, tuple, numpy.ndarray")
        def fdel(self):
            """Delete analysis variable type array"""
            del self._col_analysis_type
        return locals()
    col_analysis_type = property(**col_analysis_type())

    def col_analysis_effect():
        doc = "Analysis variable effect type {'response','fixed','random',None} array."
        def fget(self):
            """Get analysis variable effect type array"""
            if hasattr(self, "_col_analysis_effect") and isinstance(self._col_analysis_effect, dict):
                return numpy.object_(list(self._col_analysis_effect.values()))
            return None
        def fset(self, value):
            """Set analysis variable effect type array"""
            if isinstance(value, (list,tuple,numpy.ndarray)):       # if array_like
                check_len(value, "col_analysis_effect", self.ncol)  # check input length
                names = self.col_name                               # get column names
                grps = value                                        # get column grps
                value = dict(zip(names,grps))                       # construct dict of grps
            if value is None:                                       # if is None
                self._col_analysis_effect = value                   # set to None
            elif isinstance(value, dict):                           # if is dict
                check_keys_in_dict(value, "col_analysis_effect", *self._data.keys())
                check_values_in_dict_all_type(value, "col_analysis_effect", (str,type(None)))
                options = ['response', 'fixed', 'random', None]
                if any(e not in options for e in value.values()):
                    raise ValueError("unsupported value: supported values are {0}".format(options))
                self._col_analysis_effect = {k: value[k] for k in self._data.keys()}
            else:
                raise TypeError("unsupported type: supported types are dict, list, tuple, numpy.ndarray")
        def fdel(self):
            """Delete analysis variable effect type array"""
            del self._col_analysis_effect
        return locals()
    col_analysis_effect = property(**col_analysis_effect())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def col_data(self, index = None, name = None, grp = None,
                 analysis_type = None, analysis_effect = None, dtype = None,
                 return_index = False, return_name = False, return_grp = False,
                 return_analysis_type = None, return_analysis_effect = None,
                 return_dtype = False, **kwargs):
        """
        Get a column's (or columns') data from the dataframe.

        Parameters
        ----------
        index : int, None
            Integer index of the column to get.
        name : str, None
            Name of the column to get.
        grp : str, None
            Group of the column to get.
        analysis_type : str, None
            Analysis type of the column to get.
        analysis_effect : str, None
            Analysis effect type of the column to get.
        dtype : str, numpy.dtype, None
            Data type of the column to get.
        return_index : boolean, default = False
            Whether to return the column index along with the column data.
        return_name : boolean, default = False
            Whether to return the column name along with the column data.
        return_analysis_type : boolean, default = False
            Whether to return the analysis type of the column data.
        return_analysis_effect : boolean, default = False
            Whether to return the analysis effect type of the column data.
        return_grp : boolean, default = False
            Whether to return the column type along with the column data.
        return_dtype : boolean, default = False
            Whether to return the column dtype along with the column data.
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
        if grp is not None:
            mask = mask & (self.col_grp == grp)
        if analysis_type is not None:
            mask = mask & (self.col_analysis_type == analysis_type)
        if analysis_effect is not None:
            mask = mask & (self.col_analysis_effect == analysis_effect)
        if dtype is not None:
            mask = mask & (self.col_dtype == dtype)
        if mask is True:
            mask = False

        # get index of column equal to the provided name
        cix = numpy.flatnonzero(mask)

        # process errors if key is not found or multiple keys are found
        if numpy.sum(mask) == 0:
            err =  "DictDataFrame does not contain a column with matches for: \n"
            err += "    name = {0}".format(name)
            err += "    grp = {0}".format(grp)
            err += "    dtype = {0}".format(dtype)
            raise KeyError(err)

        # extract column arrays as a numpy.ndarray
        out_arr = [e for e,m in zip(self.data.values(),mask) if m]

        # construct extra output list
        out_extra = []

        # add values to extra output
        if return_index:
            out_extra.append(cix)
        if return_name:
            out_extra.append(self.col_name[cix])
        if return_grp:
            out_extra.append(self.col_grp[cix])
        if return_analysis_type:
            out_extra.append(self.analysis_type[cix])
        if return_analysis_effect:
            out_extra.append(self.analysis_effect[cix])
        if return_dtype:
            out_extra.append(self.col_dtype[cix])

        # construct output object (numpy.ndarray or tuple)
        out = out_arr if len(out_extra) == 0 else (out_arr, *out_extra)

        return out



################################################################################
################################## Utilities ###################################
################################################################################
def is_DictPhenotypeDataFrame(v):
    """
    Determine whether an object is a DictPhenotypeDataFrame.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DictPhenotypeDataFrame object instance.
    """
    return isinstance(v, DictPhenotypeDataFrame)

def check_is_DictPhenotypeDataFrame(v, vname):
    """
    Check if object is of type DictPhenotypeDataFrame. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DictPhenotypeDataFrame):
        raise TypeError("variable '{0}' must be a DictPhenotypeDataFrame".format(vname))

def cond_check_is_DictPhenotypeDataFrame(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type DictPhenotypeDataFrame. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a DictPhenotypeDataFrame.
    """
    if cond(v):
        check_is_DictPhenotypeDataFrame(v, vname)
