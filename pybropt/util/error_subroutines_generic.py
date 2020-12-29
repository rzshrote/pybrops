# generic subroutines for other error subroutines

################################################################################
########################## isinstance check functions ##########################
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

################################################################################
#################### conditional isinstance check functions ####################
################################################################################
def generic_cond_check_isinstance(v, vname, vtype, cond=(lambda s: s is not None)):
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
    if cond(v) and (not isinstance(v, vtype)):
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

################################################################################
########################### attribute check functions ##########################
################################################################################
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
    vhasattr = all(hasattr(v, e) for e in vattr) if istuple else hasattr(v, vattr)
    if not vhasattr:
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
        raise TypeError("variable '{0}' must have attribute {1}".format(vname, attrname))

################################################################################
#################### conditional attribute check functions #####################
################################################################################
def generic_cond_check_hasattr(v, vname, vattr, cond=(lambda s: s is not None)):
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
    vhasattr = all(hasattr(v, e) for e in vattr) if istuple else hasattr(v, vattr)
    if cond(v) and (not vhasattr):
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
        raise TypeError("variable '{0}' must have attribute {1}".format(vname, attrname))

################################################################################
############################## is check functions ##############################
################################################################################
def generic_check_is(v, vname, w, wname):
    if v is not w:
        raise ValueError("variable '{0}' is not '{1}'".format(vname, wname))

def generic_check_is_not(v, vname, w, wname):
    if v is w:
        raise ValueError("variable '{0}' is not '{1}'".format(vname, wname))

################################################################################
######################## conditional is check functions ########################
################################################################################
def generic_cond_check_is(v, vname, w, wname, cond=(lambda s: s is not None)):
    if cond(v) and (v is not w):
        raise ValueError("variable '{0}' is not '{1}'".format(vname, wname))

def generic_cond_check_is_not(v, vname, w, wname, cond=(lambda s: s is not None)):
    if cond(v) and (v is w):
        raise ValueError("variable '{0}' is not '{1}'".format(vname, wname))
