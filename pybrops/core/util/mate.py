"""
Module containing mating utilities.
"""

from typing import Union
import numpy
from numpy.random import Generator, RandomState

__all__ = [
    "dense_meiosis", 
    "dense_dh", 
    "dense_cross",
]

def dense_meiosis(
        geno: numpy.ndarray, 
        sel: numpy.ndarray, 
        xoprob: numpy.ndarray, 
        rng: Union[Generator,RandomState]
    ) -> numpy.ndarray:
    """
    Perform meiosis on matrix inputs. Assumes diploidy.

    Parameters
    ----------
    geno : numpy.ndarray
        A phased genome matrix of shape ``(nphase,ntaxa,nvrnt)``.

        Where:

        - ``nphase`` is the number of chromosome phases.
        - ``ntaxa`` is the number of taxa.
        - ``nvrnt`` is the number of variants (markers).
    
    sel : numpy.ndarray
        A selection configuration array of shape ``(nsel,)`` containing indices 
        of individuals for which to perform meiosis.

        Where:

        - ``nsel`` is the number of individuals for which to perform meiosis.

    xoprob : numpy.ndarray
        A crossover probability array of shape ``(nvrnt,)`` containing 
        sequential probabilities of recombination.

    rng : numpy.random.Generator, numpy.random.RandomState
        A random number generator instance.

    Returns
    -------
    gamete : numpy.ndarray
        A gamete matrix of shape ``(nsel,nvrnt)``.
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

def dense_dh(
        geno: numpy.ndarray, 
        sel: numpy.ndarray, 
        xoprob: numpy.ndarray, 
        rng: Union[Generator,RandomState]
    ) -> numpy.ndarray:
    """
    Perform doubled haploid production on matrix inputs. Assumes diploidy.

    Parameters
    ----------
    geno : numpy.ndarray
        A phased genome matrix of shape ``(nphase,ntaxa,nvrnt)``.

        Where:

        - ``nphase`` is the number of chromosome phases.
        - ``ntaxa`` is the number of taxa.
        - ``nvrnt`` is the number of variants (markers).
    
    sel : numpy.ndarray
        A selection configuration array of shape ``(nsel,)`` containing indices 
        of individuals for which to perform meiosis.

        Where:

        - ``nsel`` is the number of individuals for which to perform meiosis.

    xoprob : numpy.ndarray
        A crossover probability array of shape ``(nvrnt,)`` containing 
        sequential probabilities of recombination.

    rng : numpy.random.Generator, numpy.random.RandomState
        A random number generator instance.

    Returns
    -------
    progeny : numpy.ndarray
        A phased genome matrix of progenies of shape ``(nphase,ntaxa,nvrnt)``.
    """
    # generate gametes
    gamete = dense_meiosis(geno, sel, xoprob, rng)

    # generate offspring genotypes by stacking matrices to make 3d matrix
    progeny = numpy.stack([gamete, gamete])

    return progeny

def dense_cross(
        fgeno: numpy.ndarray, 
        mgeno: numpy.ndarray, 
        fsel: numpy.ndarray, 
        msel: numpy.ndarray, 
        xoprob: numpy.ndarray, 
        rng: Union[Generator,RandomState]
    ) -> numpy.ndarray:
    """
    Perform mating on matrix inputs.

    Parameters
    ----------
    fgeno : numpy.ndarray
        A female phased genome matrix of shape ``(nphase,ntaxa,nvrnt)``.

        Where:

        - ``nphase`` is the number of chromosome phases.
        - ``ntaxa`` is the number of taxa.
        - ``nvrnt`` is the number of variants (markers).
    
    mgeno : numpy.ndarray
        A male phased genome matrix of shape ``(nphase,ntaxa,nvrnt)``.

        Where:

        - ``nphase`` is the number of chromosome phases.
        - ``ntaxa`` is the number of taxa.
        - ``nvrnt`` is the number of variants (markers).
    
    fsel : numpy.ndarray
        A female selection configuration array of shape ``(nsel,)`` containing 
        indices of individuals for which to perform meiosis.

        Where:

        - ``nsel`` is the number of individuals for which to perform meiosis.

    msel : numpy.ndarray
        A selection configuration array of shape ``(nsel,)`` containing indices 
        of individuals for which to perform meiosis.

        Where:

        - ``nsel`` is the number of individuals for which to perform meiosis.

    xoprob : numpy.ndarray
        A crossover probability array of shape ``(nvrnt,)`` containing 
        sequential probabilities of recombination.

    rng : numpy.random.Generator, numpy.random.RandomState
        A random number generator instance.

    Returns
    -------
    progeny : numpy.ndarray
        Genotype matrix of progenies.
    """
    # generate gametes
    fgamete = dense_meiosis(fgeno, fsel, xoprob, rng)
    mgamete = dense_meiosis(mgeno, msel, xoprob, rng)

    # generate offspring genotypes by stacking matrices to make 3d matrix
    progeny = numpy.stack([fgamete, mgamete])

    return progeny
