# generic subroutines for other error subroutines

################################################################################
############################### check functions ################################
################################################################################
def generic_check_isinstance(v, vname, vtype):
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

def generic_check_hasattr(v, vname, vattr):
    """
    Generic check has attribute function.

    Parameters
    ----------
    v : object
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vattr : str, tuple
        type : Attribute string associated with the object variable.
        tuple : Tuple of attribute strings (logical and) for the object variable.
    """
    istuple = isinstance(vattr, tuple)
    logic = all(hasattr(v, e) for e in vattr) if istuple else hasattr(v, vattr)
    if not logic:
        attrname = None
        if istuple:
            attrname = ""
            l = len(vtype)
            for e, i in enumerate(vtype):
                attrname += e
                if i < l - 2:
                    attrname += ", "
                elif i < l - 1:
                    attrname += ", and "
        else:
            attrname = vattr
        raise AttributeError("variable '{0}' must have attribute {1}".format(vname, attrname))

def generic_check_is(v, vname, w, wname):
    if v is not w:
        raise ValueError("variable '{0}' is not '{1}'".format(vname, wname))

def generic_check_is_not(v, vname, w, wname):
    if v is w:
        raise ValueError("variable '{0}' is not '{1}'".format(vname, wname))

def generic_check_len(v, vname, n):
    if len(v) != n:
        raise ValueError("the length of '{0}' is not equal to {1}".format(vname,n))

def generic_check_len_eq(v, vname, w, wname):
    if len(v) != len(w):
        raise ValueError("the lengths of '{0}' and '{1}' are not equivalent".format(vname, wname))

def generic_check_dict_keys(v, vname, w, wname):
    if any(e not in v for e in w):
        raise ValueError("dict '{0}' must have keys: {1}".format(vname, wname))

def generic_check_gteq(v, vname, w):
    if v < w:
        raise ValueError("variable '{0}' is not greater than or equal to {1}".format(vname, w))

################################################################################
######################### conditional check functions ##########################
################################################################################
def generic_default_cond(s):
    return s is not None

def generic_cond_check_isinstance(v, vname, vtype, cond = generic_default_cond):
    """
    Generic conditional check type function.

    Parameters
    ----------
    v : object
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vtype : type, tuple
        Type associated with the object variable or tuple of acceptable types
        the object variable can be.
    cond : callable
        Function to determine whether or not to check variable type.
    """
    if cond(v):
        generic_check_isinstance(v, vname, vtype)

def generic_cond_check_hasattr(v, vname, vattr, cond = generic_default_cond):
    """
    Generic conditional check has attribute function.

    Parameters
    ----------
    v : object
        Reference to object variable.
    vname : str
        Name associated with the object variable.
    vattr : str, tuple
        type : Attribute string associated with the object variable.
        tuple : Tuple of attribute strings (logical and) for the object variable.
    """
    if cond(v):
        generic_check_hasattr(v, vname, vattr)

def generic_cond_check_is(v, vname, w, wname, cond = generic_default_cond):
    if cond(v):
        generic_check_is(v, vname, w, wname)

def generic_cond_check_is_not(v, vname, w, wname, cond = generic_default_cond):
    if cond(v):
        generic_check_is_not(v, vname, w, wname)

def generic_cond_check_len(v, vname, n, cond = generic_default_cond):
    if cond(v):
        generic_check_len(v, vname, n)

def generic_cond_check_len_eq(v, vname, w, wname, cond = generic_default_cond):
    if cond(v):
        generic_check_len_eq(v, vname, w, wname)

def generic_cond_check_dict_keys(v, vname, w, wname, cond = generic_default_cond):
    if cond(v):
        generic_check_dict_keys(v, vname, w, wname)

def generic_cond_check_gteq(v, vname, w, cond = generic_default_cond):
    if cond(v):
        generic_check_gteq(v, vname, w)
