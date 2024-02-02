"""
Module containing mating utilities.
"""

from typing import Union
import numpy

__all__ = [
    "mat_meiosis", 
    "mat_dh", 
    "mat_mate",
]

def mat_meiosis(
        geno: numpy.ndarray, 
        sel: numpy.ndarray, 
        xoprob: numpy.ndarray, 
        rng: Union[numpy.random.Generator,numpy.random.RandomState]
    ) -> numpy.ndarray:
    """
    Perform meiosis on matrix inputs.

    Parameters
    ----------
    geno : numpy.ndarray
        Genotype matrix.
    sel : numpy.ndarray
        Selection configuration array.
    xoprob : numpy.ndarray
        Crossover porbability array.
    rng : numpy.random.Generator, numpy.random.RandomState
        Random number generator instance

    Returns
    -------
    gamete : numpy.ndarray
        Genotype matrix of gametes.
    """
    # calculate shape of the retured gamete matrix
    # number of rows is the number of elements in sel
    gshape = (len(sel), len(xoprob))

    # generate random numbers to determine crossover points:
    # generate in interval [0,1)
    rnd = rng.uniform(0, 1, gshape)

    # allocate gamete array
    gamete = numpy.empty(gshape, dtype = geno.dtype)

    for i,s in enumerate(sel):
        # calculate locations where crossover occurs
        xoix = numpy.flatnonzero(rnd[i] < xoprob)

        # get starting phase
        phase = 0

        # starting index for copy
        stix = 0

        for spix in xoix:
            # copy from start index (inclusive) to stop index (exclusive)
            gamete[i,stix:spix] = geno[phase,s,stix:spix]

            # move the start index in the next iteration to the current stop index
            stix = spix

            # alternate between phases
            phase = 1 - phase

        # finally copy the last remaining chromosome segments
        gamete[i,stix:] = geno[phase,s,stix:]

    return gamete

def mat_dh(
        geno: numpy.ndarray, 
        sel: numpy.ndarray, 
        xoprob: numpy.ndarray, 
        rng: Union[numpy.random.Generator,numpy.random.RandomState]
    ) -> numpy.ndarray:
    """
    Perform doubled haploid production on matrix inputs.

    Parameters
    ----------
    geno : numpy.ndarray
        Genotype matrix.
    sel : numpy.ndarray
        Selection configuration array.
    xoprob : numpy.ndarray
        Crossover porbability array.
    rng : numpy.random.Generator, numpy.random.RandomState
        Random number generator instance

    Returns
    -------
    progeny : numpy.ndarray
        Genotype matrix of progenies.
    """
    # generate gametes
    gamete = mat_meiosis(geno, sel, xoprob, rng)

    # generate offspring genotypes by stacking matrices to make 3d matrix
    progeny = numpy.stack([gamete, gamete])

    return progeny

def mat_mate(
        fgeno: numpy.ndarray, 
        mgeno: numpy.ndarray, 
        fsel: numpy.ndarray, 
        msel: numpy.ndarray, 
        xoprob: numpy.ndarray, 
        rng: Union[numpy.random.Generator,numpy.random.RandomState]
    ) -> numpy.ndarray:
    """
    Perform mating on matrix inputs.

    Parameters
    ----------
    fgeno : numpy.ndarray
        Female genotype matrix.
    mgeno : numpy.ndarray
        Male genotype matrix.
    fsel : numpy.ndarray
        Female selection configuration array.
    msel : numpy.ndarray
        Male selection configuration array.
    xoprob : numpy.ndarray
        Crossover porbability array.
    rng : numpy.random.Generator, numpy.random.RandomState
        Random number generator instance

    Returns
    -------
    progeny : numpy.ndarray
        Genotype matrix of progenies.
    """
    # generate gametes
    fgamete = mat_meiosis(fgeno, fsel, xoprob, rng)
    mgamete = mat_meiosis(mgeno, msel, xoprob, rng)

    # generate offspring genotypes by stacking matrices to make 3d matrix
    progeny = numpy.stack([fgamete, mgamete])

    return progeny
