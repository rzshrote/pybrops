import numbers

from . import generic_check_isinstance

from . import generic_default_cond
from . import generic_cond_check_isinstance

################################################################################
###################### basic inheritance check functions #######################
################################################################################
def check_inherits(obj, objname, objtype):
    """
    Generic check of inheritance using method resolution order metadata.

    Parameters
    ----------
    obj : object
        Python object to check.
    objname : str
        Name associated with the object variable.
    objtype : type
        Object type from which 'obj' must inherit.
    """
    if isinstance(objtype, type) and (objname not in obj.__mro__):
        raise TypeError("variable '{0}' must inherit from '{1}'".format(objname, objtype))
    elif isinstance(objtype, tuple) and (all(e not in obj.__mro__ for e in objtype)):
        raise TypeError("variable '{0}' must inherit from one of '{1}'".format(objname, objtype))
    else:
        raise TypeError("'objtype' must be of type 'type' or 'tuple'")

################################################################################
################## basic check functions for basic data types ##################
################################################################################
def check_isinstance(v, vname, vtype):
    """
    Generic check type function.

    Parameters
    ----------
    v : object
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vtype : type, tuple
        type : Type associated with the object variable.
        tuple : Tuple of acceptable types (logical or) for the object variable.
    """
    if not isinstance(v, vtype):
        tname = None
        if isinstance(vtype, tuple):
            tname = ""
            l = len(vtype)
            for e, i in enumerate(vtype):
                tname += e.__name__
                if i < l - 2:
                    tname += ", "
                elif i < l - 1:
                    tname += ", or "
        else:
            tname = vtype.__name__
        raise TypeError("variable '{0}' must of type {1}".format(vname, tname))

def check_is_bool(v, vname):
    generic_check_isinstance(v, vname, bool)

# def check_is_bytearray(v, vname):
#     generic_check_isinstance(v, vname, bytearray)

# def check_is_bytes(v, vname):
#     generic_check_isinstance(v, vname, bytes)

# def check_is_complex(v, vname):
#     generic_check_isinstance(v, vname, complex)

def check_is_dict(v, vname):
    generic_check_isinstance(v, vname, dict)

def check_is_float(v, vname):
    generic_check_isinstance(v, vname, float)

# def check_is_frozenset(v, vname):
#     generic_check_isinstance(v, vname, frozenset)

def check_is_int(v, vname):
    generic_check_isinstance(v, vname, int)

def check_is_list(v, vname):
    generic_check_isinstance(v, vname, list)

# def check_is_memoryview(v, vname):
#     generic_check_isinstance(v, vname, memoryview)

def check_is_range(v, vname):
    generic_check_isinstance(v, vname, range)

def check_is_set(v, vname):
    generic_check_isinstance(v, vname, set)

def check_is_str(v, vname):
    generic_check_isinstance(v, vname, str)

def check_is_tuple(v, vname):
    generic_check_isinstance(v, vname, tuple)

def check_is_type(obj, objname):
    if not isinstance(obj, type):
        raise TypeError("variable '{0}' must be of type 'type'".format(objname))

def check_is_Number(v, vname):
    generic_check_isinstance(v, vname, numbers.Number)

def check_is_Integral(v, vname):
    generic_check_isinstance(v, vname, numbers.Integral)

################################################################################
################ compound check functions for basic data types #################
################################################################################
def check_is_array_like(v, vname):
    alattr = ("__len__","__iter__","__getitem__")
    for a in alattr:
        if not hasattr(v, a):
            raise TypeError("'{0}' must have attribute '{1}' to be array_like".format(vname,a))

def check_is_str_or_iterable(v, vname):
    if not (isinstance(s, str) or hasattr(s, "__iter__")):
        raise TypeError("'{0}' must be of type str or have attribute __iter__.".format(vname))

def check_is_list_or_tuple(v, vname):
    generic_check_isinstance(v, vname, (list,tuple))


################################################################################
######################### conditional check functions ##########################
################################################################################
def cond_check_is_bool(v, vname):
    generic_cond_check_isinstance(v, vname, bool)

# def cond_check_is_bytearray(v, vname):
#     generic_cond_check_isinstance(v, vname, bytearray)

# def cond_check_is_bytes(v, vname):
#     generic_cond_check_isinstance(v, vname, bytes)

# def cond_check_is_complex(v, vname):
#     generic_cond_check_isinstance(v, vname, complex)

def cond_check_is_dict(v, vname):
    generic_cond_check_isinstance(v, vname, dict)

def cond_check_is_float(v, vname):
    generic_cond_check_isinstance(v, vname, float)

# def cond_check_is_frozenset(v, vname):
#     generic_cond_check_isinstance(v, vname, frozenset)

def cond_check_is_int(v, vname):
    generic_cond_check_isinstance(v, vname, int)

def cond_check_is_list(v, vname):
    generic_cond_check_isinstance(v, vname, list)

# def cond_check_is_memoryview(v, vname):
#     generic_cond_check_isinstance(v, vname, memoryview)

def cond_check_is_range(v, vname):
    generic_cond_check_isinstance(v, vname, range)

def cond_check_is_set(v, vname):
    generic_cond_check_isinstance(v, vname, set)

def cond_check_is_str(v, vname):
    generic_cond_check_isinstance(v, vname, str)

def cond_check_is_tuple(v, vname):
    generic_cond_check_isinstance(v, vname, tuple)


################################################################################
########## conditional compound check functions for basic data types ###########
################################################################################

def cond_check_is_str_or_iterable(v, vname, cond = generic_default_cond):
    if cond(v):
        check_is_str_or_iterable(v, vname)