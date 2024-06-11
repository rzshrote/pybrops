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
    "twowayix_symab_anyab_len",
    "twowayix_symab_selfab_len",
    "twowayix_symab_uniqab_len",
    "twowayix_asymab_anyab_len",
    "twowayix_asymab_selfab_len",
    "twowayix_asymab_uniqab_len",
    "twowayix_len",
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

######################## Generators ########################

def twowayix_symab_anyab(ntaxa: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow selfing, allow outcrossing

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if ntaxa <= 0:
        yield from ()
    for female in range(ntaxa):
        for male in range(female,ntaxa):
            yield (female,male)

def twowayix_symab_selfab(ntaxa: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow selfing, disallow outcrossing

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if ntaxa <= 0:
        yield from ()
    for female in range(ntaxa):
        yield (female,female)

def twowayix_symab_uniqab(ntaxa: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Disallow selfing, allow outcrossing

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if ntaxa <= 1:
        yield from ()
    for female in range(ntaxa):
        for male in range(female+1,ntaxa):
            yield (female,male)

def twowayix_asymab_anyab(ntaxa: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow selfing, allow outcrossing

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if ntaxa <= 0:
        yield from ()
    for female in range(ntaxa):
        for male in range(ntaxa):
            yield (female,male)

def twowayix_asymab_selfab(ntaxa: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow selfing, disallow outcrossing

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if ntaxa <= 0:
        yield from ()
    for female in range(ntaxa):
        yield (female,female)

def twowayix_asymab_uniqab(ntaxa: int) -> Generator:
    """
    Generate indices for two-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Disallow selfing, allow outcrossing

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (female,male) == (A x B).
    """
    if ntaxa <= 0:
        yield from ()
    for female in range(ntaxa):
        for male in range(female):
            yield (female,male)
        for male in range(female+1,ntaxa):
            yield (female,male)

def twowayix(ntaxa: int, symab: bool = True, mateab: str = "uniq") -> Generator:
    """
    Generate two-way cross indices.

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

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
            yield from twowayix_symab_anyab(ntaxa)
        elif mateab == "self":
            yield from twowayix_symab_selfab(ntaxa)
        elif mateab == "uniq":
            yield from twowayix_symab_uniqab(ntaxa)
        else:
            raise ValueError
    else:
        if mateab == "any":
            yield from twowayix_asymab_anyab(ntaxa)
        elif mateab == "self":
            yield from twowayix_asymab_selfab(ntaxa)
        elif mateab == "uniq":
            yield from twowayix_asymab_uniqab(ntaxa)
        else:
            raise ValueError

################### Length of Generators ###################

def twowayix_symab_anyab_len(ntaxa: int) -> int:
    """
    Calculate the length of a ``twowayix_symab_anyab`` generator.

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Returns
    -------
    out : int
        The length of the ``twowayix_symab_anyab`` generator.
    """
    if ntaxa <= 0:
        return 0
    out = (ntaxa * (ntaxa+1))//2
    return out

def twowayix_symab_selfab_len(ntaxa: int) -> int:
    """
    Calculate the length of a ``twowayix_symab_selfab`` generator.

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Returns
    -------
    out : int
        The length of the ``twowayix_symab_selfab`` generator.
    """
    if ntaxa <= 0:
        return 0
    out = ntaxa
    return out

def twowayix_symab_uniqab_len(ntaxa: int) -> int:
    """
    Calculate the length of a ``twowayix_symab_uniqab`` generator.

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Returns
    -------
    out : int
        The length of the ``twowayix_symab_uniqab`` generator.
    """
    if ntaxa <= 1:
        return 0
    out = (ntaxa * (ntaxa-1))//2
    return out

def twowayix_asymab_anyab_len(ntaxa: int) -> int:
    """
    Calculate the length of a ``twowayix_asymab_anyab`` generator.

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Returns
    -------
    out : int
        The length of the ``twowayix_asymab_anyab`` generator.
    """
    if ntaxa <= 0:
        return 0
    out = ntaxa * ntaxa
    return out

def twowayix_asymab_selfab_len(ntaxa: int) -> int:
    """
    Calculate the length of a ``twowayix_asymab_selfab`` generator.

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Returns
    -------
    out : int
        The length of the ``twowayix_asymab_selfab`` generator.
    """
    if ntaxa <= 0:
        return 0
    out = ntaxa
    return out

def twowayix_asymab_uniqab_len(ntaxa: int) -> int:
    """
    Calculate the length of a ``twowayix_asymab_uniqab`` generator.

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Returns
    -------
    out : int
        The length of the ``twowayix_asymab_uniqab`` generator.
    """
    if ntaxa <= 0:
        return 0
    out = ntaxa * (ntaxa - 1)
    return out

def twowayix_len(ntaxa: int, symab: bool = True, mateab: str = "uniq") -> int:
    """
    Calculate the length of a ``twowayix`` generator.

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    symab : bool, default = True
        Whether a cross is symmetric or not.
        If ``True``, assume cross symmetry ( (AxB) == (BxA) ).
        If ``False``, assumer cross asymmetry ( (AxB) != (BxA) ).
    
    mateab : str, default = "uniq"
        Mating strategy. Must be one of {"any", "self", "uniq"}.
        If ``"any"``, no parental restrictions (selfing, outcrossing allowed).
        If ``"self"``, require identical parents (selfing allowed; outcrossing disallowed).
        If ``"uniq"``, require unique parents (outcrossing allowed; selfing disallowed)

    Returns
    -------
    out : int
        The length of the ``twowayix`` generator.
    """
    if symab:
        if mateab == "any":
            return twowayix_symab_anyab_len(ntaxa)
        elif mateab == "self":
            return twowayix_symab_selfab_len(ntaxa)
        elif mateab == "uniq":
            return twowayix_symab_uniqab_len(ntaxa)
        else:
            raise ValueError
    else:
        if mateab == "any":
            return twowayix_asymab_anyab_len(ntaxa)
        elif mateab == "self":
            return twowayix_asymab_selfab_len(ntaxa)
        elif mateab == "uniq":
            return twowayix_asymab_uniqab_len(ntaxa)
        else:
            raise ValueError

####################### Three-way cross index generators #######################

def threewayix_symab_anyab_anyc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for recurrent in range(ntaxa):
        for female in range(ntaxa):
            for male in range(female,ntaxa):
                yield (recurrent,female,male)

def threewayix_symab_anyab_backc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for recurrent in range(ntaxa):
        female = recurrent
        for male in range(female,ntaxa):
            yield (recurrent,female,male)
    # for i in range(ntaxa):
    #     for j in range(i,ntaxa):
    #         yield (i,i,j)
    #         yield (j,i,j)

def threewayix_symab_anyab_uniqc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(i,ntaxa):
            for k in range(i):
                yield (k,i,j)
            for k in range(i+1,j):
                yield (k,i,j)
            for k in range(k+1,ntaxa):
                yield (k,i,j)

def threewayix_symab_selfab_anyc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(ntaxa):
            yield (j,i,i)

def threewayix_symab_selfab_backc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        yield (i,i,i)

def threewayix_symab_selfab_uniqc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(i):
            yield (j,i,i)
        for j in range(i+1,ntaxa):
            yield (j,i,i)

def threewayix_symab_uniqab_anyc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(i+1,ntaxa):
            for k in range(ntaxa):
                yield (k,i,j)

def threewayix_symab_uniqab_backc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(i+1,ntaxa):
            yield (i,i,j)
            yield (j,i,j)

def threewayix_symab_uniqab_uniqc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume symmetric mate pairings for (AxB) cross: (A x B) == (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(i+1,ntaxa):
            for k in range(i):
                yield (k,i,j)
            for k in range(i+1,j):
                yield (k,i,j)
            for k in range(k+1,ntaxa):
                yield (k,i,j)

def threewayix_asymab_anyab_anyc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(ntaxa):
            for k in range(ntaxa):
                yield (k,i,j)

def threewayix_asymab_anyab_backc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(ntaxa):
            yield (i,i,j)
            yield (j,i,j)

def threewayix_asymab_anyab_uniqc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(ntaxa):
            for k in range(i):
                yield (k,i,j)
            for k in range(i+1,j):
                yield (k,i,j)
            for k in range(j+1,ntaxa):
                yield (k,i,j)

def threewayix_asymab_selfab_anyc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(ntaxa):
            yield (j,i,i)

def threewayix_asymab_selfab_backc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        yield (i,i,i)

def threewayix_asymab_selfab_uniqc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Allow (AxB) selfing, disallow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(i):
            yield (j,i,i)
        for j in range(i+1,ntaxa):
            yield (j,i,i)

def threewayix_asymab_uniqab_anyc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Any parent can be the recurrent parent: C == (A, B, any).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(i):
            for k in range(ntaxa):
                yield (k,i,j)
        for j in range(i+1,ntaxa):
            for k in range(ntaxa):
                yield (k,i,j)

def threewayix_asymab_uniqab_backc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit backcrossing for the recurrent parent: C == (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(i):
            yield (i,i,j)
            yield (j,i,j)
        for j in range(i+1,ntaxa):
            yield (i,i,j)
            yield (j,i,j)

def threewayix_asymab_uniqab_uniqc(ntaxa: int) -> Generator:
    """
    Generate indices for three-way parent crosses with the following constraints:

    1) Assume asymmetric mate pairings for (AxB) cross: (A x B) != (B x A)
    2) Disallow (AxB) selfing, allow (AxB) outcrossing
    3) Only permit outcrossing for the recurrent parent: C != (A or B).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices. Indices are (recurrent,female,male) == (C x (A x B)).
    """
    for i in range(ntaxa):
        for j in range(i):
            for k in range(i):
                yield (k,i,j)
            for k in range(i+1,j):
                yield (k,i,j)
            for k in range(k+1,ntaxa):
                yield (k,i,j)
        for j in range(i+1,ntaxa):
            for k in range(i):
                yield (k,i,j)
            for k in range(i+1,j):
                yield (k,i,j)
            for k in range(k+1,ntaxa):
                yield (k,i,j)

def threewayix(ntaxa: int, symab: bool = True, mateab: str = "uniq", matec: str = "uniq") -> Generator:
    """
    Generate three-way cross indices.

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

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
                yield from threewayix_symab_anyab_anyc(ntaxa)
            elif matec == "back":
                yield from threewayix_symab_anyab_backc(ntaxa)
            elif matec == "uniq":
                yield from threewayix_symab_anyab_uniqc(ntaxa)
            else:
                raise ValueError
        elif mateab == "self":
            if matec == "any":
                yield from threewayix_symab_selfab_anyc(ntaxa)
            elif matec == "back":
                yield from threewayix_symab_selfab_backc(ntaxa)
            elif matec == "uniq":
                yield from threewayix_symab_selfab_uniqc(ntaxa)
            else:
                raise ValueError
        elif mateab == "uniq":
            if matec == "any":
                yield from threewayix_symab_uniqab_anyc(ntaxa)
            elif matec == "back":
                yield from threewayix_symab_uniqab_backc(ntaxa)
            elif matec == "uniq":
                yield from threewayix_symab_uniqab_uniqc(ntaxa)
            else:
                raise ValueError
        else:
            raise ValueError
    else:
        if mateab == "any":
            if matec == "any":
                yield from threewayix_asymab_anyab_anyc(ntaxa)
            elif matec == "back":
                yield from threewayix_asymab_anyab_backc(ntaxa)
            elif matec == "uniq":
                yield from threewayix_asymab_anyab_uniqc(ntaxa)
            else:
                raise ValueError
        elif mateab == "self":
            if matec == "any":
                yield from threewayix_asymab_selfab_anyc(ntaxa)
            elif matec == "back":
                yield from threewayix_asymab_selfab_backc(ntaxa)
            elif matec == "uniq":
                yield from threewayix_asymab_selfab_uniqc(ntaxa)
            else:
                raise ValueError
        elif mateab == "uniq":
            if matec == "any":
                yield from threewayix_asymab_uniqab_anyc(ntaxa)
            elif matec == "back":
                yield from threewayix_asymab_uniqab_backc(ntaxa)
            elif matec == "uniq":
                yield from threewayix_asymab_uniqab_uniqc(ntaxa)
            else:
                raise ValueError
        else:
            raise ValueError

####################### Four-way cross index generators ########################

def fourwayix_perm_self(ntaxa: int) -> Generator:
    """
    Generate all four-way parent permutations (order matters).

    Parameters
    ----------
    ntaxa : int
        Number of taxa eligible to serve as parents.

    Yields
    ------
    out : tuple
        Tuple of indices.
    """
    for i in range(ntaxa):
        for j in range(ntaxa):
            for k in range(ntaxa):
                for l in range(ntaxa):
                    yield (i,j,k,l)
