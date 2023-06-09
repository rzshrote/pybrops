"""
Module containing generic subroutines for other error subroutines.
"""

################################################################################
############################### check functions ################################
################################################################################


def generic_check_isinstance(v: object, vname: str, vtype: type) -> None:
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
        raise TypeError("variable '{0}' must be of type {1}".format(vname, tname))
