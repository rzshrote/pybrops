"""
Module containing utility functions to generate cross indices.
"""

__all__ = [
    # two-way crosses
    "twowayix_symab_anyab",
    "twowayix_symab_selfab",
    "twowayix_symab_uniqab",
    "twowayix_asymab_anyab",
    "twowayix_asymab_selfab",
    "twowayix_asymab_uniqab",
    "twowayix",
    # three-way crosses
    "threewayix_symab_anyab_anyc",
    "threewayix_symab_anyab_backc",
    "threewayix_symab_anyab_uniqc",
    "threewayix_symab_selfab_anyc",
    "threewayix_symab_selfab_backc",
    "threewayix_symab_selfab_uniqc",
    "threewayix_symab_uniqab_anyc",
    "threewayix_symab_uniqab_backc",
    "threewayix_symab_uniqab_uniqc",
    "threewayix_asymab_anyab_anyc",
    "threewayix_asymab_anyab_backc",
    "threewayix_asymab_anyab_uniqc",
    "threewayix_asymab_selfab_anyc",
    "threewayix_asymab_selfab_backc",
    "threewayix_asymab_selfab_uniqc",
    "threewayix_asymab_uniqab_anyc",
    "threewayix_asymab_uniqab_backc",
    "threewayix_asymab_uniqab_uniqc",
    "threewayix",
    # four-way crosses
    "fourwayix_perm_self",
]

from typing import Generator

######################## Two-way cross index generators ########################

def twowayix_symab_anyab(nparent: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow selfing, allow outcrossing

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if nparent <= 0:
        yield from ()
    for i in range(nparent):
        for j in range(i,nparent):
            yield (i,j)

def twowayix_symab_selfab(nparent: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow selfing, disallow outcrossing

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if nparent <= 0:
        yield from ()
    for i in range(nparent):
        yield (i,i)

def twowayix_symab_uniqab(nparent: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Disallow selfing, allow outcrossing

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if nparent <= 0:
        yield from ()
    for i in range(nparent):
        for j in range(i+1,nparent):
            yield (i,j)

def twowayix_asymab_anyab(nparent: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow selfing, allow outcrossing

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if nparent <= 0:
        yield from ()
    for i in range(nparent):
        for j in range(nparent):
            yield (i,j)

def twowayix_asymab_selfab(nparent: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow selfing, disallow outcrossing

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if nparent <= 0:
        yield from ()
    for i in range(nparent):
        yield (i,i)

def twowayix_asymab_uniqab(nparent: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Disallow selfing, allow outcrossing

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if nparent <= 0:
        yield from ()
    for i in range(nparent):
        for j in range(i):
            yield (i,j)
        for j in range(i+1,nparent):
            yield (i,j)

def twowayix(nparent: int, symab: bool = True, mateab: str = "uniq") -> Generator:
    """
    Generate two-way cross indices.

    Parameters
    ----------
    nparent : int
        Number of parents.

    symab : bool, default = True
        Whether a cross is symmetric or not.
        If ``True``, assume cross symmetry ( (AxB) == (BxA) ).
        If ``False``, assumer cross asymmetry ( (AxB) != (BxA) ).
    
    mateab : str, default = "uniq"
        Mating strategy. Must be one of {"any", "self", "uniq"}.
        If ``"any"``, no parental restrictions (selfing, outcrossing allowed).
        If ``"self"``, require identical parents (selfing allowed; outcrossing disallowed).
        If ``"uniq"``, require unique parents (outcrossing allowed; selfing disallowed)

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if symab:
        if mateab == "any":
            yield from twowayix_symab_anyab(nparent)
        elif mateab == "self":
            yield from twowayix_symab_selfab(nparent)
        elif mateab == "uniq":
            yield from twowayix_symab_uniqab(nparent)
        else:
            raise ValueError
    else:
        if mateab == "any":
            yield from twowayix_asymab_anyab(nparent)
        elif mateab == "self":
            yield from twowayix_asymab_selfab(nparent)
        elif mateab == "uniq":
            yield from twowayix_asymab_uniqab(nparent)
        else:
            raise ValueError

####################### Three-way cross index generators #######################

def threewayix_symab_anyab_anyc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(i,nparent):
            for k in range(nparent):
                yield (k,i,j)

def threewayix_symab_anyab_backc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(i,nparent):
            yield (i,i,j)
            yield (j,i,j)

def threewayix_symab_anyab_uniqc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(i,nparent):
            for k in range(i):
                yield (k,i,j)
            for k in range(i+1,j):
                yield (k,i,j)
            for k in range(k+1,nparent):
                yield (k,i,j)

def threewayix_symab_selfab_anyc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(nparent):
            yield (j,i,i)

def threewayix_symab_selfab_backc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        yield (i,i,i)

def threewayix_symab_selfab_uniqc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(i):
            yield (j,i,i)
        for j in range(i+1,nparent):
            yield (j,i,i)

def threewayix_symab_uniqab_anyc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(i+1,nparent):
            for k in range(nparent):
                yield (k,i,j)

def threewayix_symab_uniqab_backc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(i+1,nparent):
            yield (i,i,j)
            yield (j,i,j)

def threewayix_symab_uniqab_uniqc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(i+1,nparent):
            for k in range(i):
                yield (k,i,j)
            for k in range(i+1,j):
                yield (k,i,j)
            for k in range(k+1,nparent):
                yield (k,i,j)

def threewayix_asymab_anyab_anyc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(nparent):
            for k in range(nparent):
                yield (k,i,j)

def threewayix_asymab_anyab_backc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(nparent):
            yield (i,i,j)
            yield (j,i,j)

def threewayix_asymab_anyab_uniqc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(nparent):
            for k in range(i):
                yield (k,i,j)
            for k in range(i+1,j):
                yield (k,i,j)
            for k in range(j+1,nparent):
                yield (k,i,j)

def threewayix_asymab_selfab_anyc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(nparent):
            yield (j,i,i)

def threewayix_asymab_selfab_backc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        yield (i,i,i)

def threewayix_asymab_selfab_uniqc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(i):
            yield (j,i,i)
        for j in range(i+1,nparent):
            yield (j,i,i)

def threewayix_asymab_uniqab_anyc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(i):
            for k in range(nparent):
                yield (k,i,j)
        for j in range(i+1,nparent):
            for k in range(nparent):
                yield (k,i,j)

def threewayix_asymab_uniqab_backc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(i):
            yield (i,i,j)
            yield (j,i,j)
        for j in range(i+1,nparent):
            yield (i,i,j)
            yield (j,i,j)

def threewayix_asymab_uniqab_uniqc(nparent: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(nparent):
        for j in range(i):
            for k in range(i):
                yield (k,i,j)
            for k in range(i+1,j):
                yield (k,i,j)
            for k in range(k+1,nparent):
                yield (k,i,j)
        for j in range(i+1,nparent):
            for k in range(i):
                yield (k,i,j)
            for k in range(i+1,j):
                yield (k,i,j)
            for k in range(k+1,nparent):
                yield (k,i,j)

def threewayix(nparent: int, symab: bool = True, mateab: str = "uniq", matec: str = "uniq") -> Generator:
    """
    Generate three-way cross indices.

    Parameters
    ----------
    nparent : int
        Number of parents.

    symab : bool, default = True
        Whether the AxB cross is symmetric or not.
        If ``True``, assume cross symmetry ( (AxB) == (BxA) ).
        If ``False``, assumer cross asymmetry ( (AxB) != (BxA) ).
    
    mateab : str, default = "uniq"
        Mating strategy for AxB cross. Must be one of {"any", "self", "uniq"}.
        If ``"any"``, no parental restrictions (selfing, outcrossing allowed).
        If ``"self"``, require identical parents (selfing allowed; outcrossing disallowed).
        If ``"uniq"``, require unique parents (outcrossing allowed; selfing disallowed)

    matec : str, default = "uniq"
        Mating strategy for cross with C. Must be one of {"any", "back", "uniq"}.
        If ``"any"``, no parental restrictions (backcrossing, outcrossing allowed).
        If ``"back"``, require backcross parents (backcrossing allowed; outcrossing disallowed).
        If ``"uniq"``, require unique parents (outcrossing allowed; backcrossing disallowed)

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    if symab:
        if mateab == "any":
            if matec == "any":
                yield from threewayix_symab_anyab_anyc(nparent)
            elif matec == "back":
                yield from threewayix_symab_anyab_backc(nparent)
            elif matec == "uniq":
                yield from threewayix_symab_anyab_uniqc(nparent)
            else:
                raise ValueError
        elif mateab == "self":
            if matec == "any":
                yield from threewayix_symab_selfab_anyc(nparent)
            elif matec == "back":
                yield from threewayix_symab_selfab_backc(nparent)
            elif matec == "uniq":
                yield from threewayix_symab_selfab_uniqc(nparent)
            else:
                raise ValueError
        elif mateab == "uniq":
            if matec == "any":
                yield from threewayix_symab_uniqab_anyc(nparent)
            elif matec == "back":
                yield from threewayix_symab_uniqab_backc(nparent)
            elif matec == "uniq":
                yield from threewayix_symab_uniqab_uniqc(nparent)
            else:
                raise ValueError
        else:
            raise ValueError
    else:
        if mateab == "any":
            if matec == "any":
                yield from threewayix_asymab_anyab_anyc(nparent)
            elif matec == "back":
                yield from threewayix_asymab_anyab_backc(nparent)
            elif matec == "uniq":
                yield from threewayix_asymab_anyab_uniqc(nparent)
            else:
                raise ValueError
        elif mateab == "self":
            if matec == "any":
                yield from threewayix_asymab_selfab_anyc(nparent)
            elif matec == "back":
                yield from threewayix_asymab_selfab_backc(nparent)
            elif matec == "uniq":
                yield from threewayix_asymab_selfab_uniqc(nparent)
            else:
                raise ValueError
        elif mateab == "uniq":
            if matec == "any":
                yield from threewayix_asymab_uniqab_anyc(nparent)
            elif matec == "back":
                yield from threewayix_asymab_uniqab_backc(nparent)
            elif matec == "uniq":
                yield from threewayix_asymab_uniqab_uniqc(nparent)
            else:
                raise ValueError
        else:
            raise ValueError

####################### Four-way cross index generators ########################

def fourwayix_perm_self(nparent: int) -> Generator:
    """
    Generate all four-way parent permutations (order matters).

    Parameters
    ----------
    nparent : int
        Number of parents.

    Yields
    ------
    out : tuple
        Tuple of indices.
    """
    for i in range(nparent):
        for j in range(nparent):
            for k in range(nparent):
                for l in range(nparent):
                    yield (i,j,k,l)
