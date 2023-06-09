"""
Module containing haplotype calculation utility functions.
"""

from numbers import Integral
import numpy

__all__ = [
    "nhaploblk_chrom", 
    "haplobin", 
    "haplobin_bounds",
    "haplomat"
]

def nhaploblk_chrom(
        nhaploblk: int, 
        genpos: numpy.ndarray, 
        chrgrp_stix: numpy.ndarray, 
        chrgrp_spix: numpy.ndarray
    ) -> numpy.ndarray:
    """
    Given a total number of haplotype blocks to assign across the genome,
    determine the number of haplotype blocks to give to each chromosome.

    Parameters
    ----------
    nhaploblk : int
        Total number of haplotype blocks to assign across the genome.
    genpos : numpy.ndarray
        Array of shape ``(p,)`` containing genetic positions of markers across
        the genome.

        Where:

        - ``p`` is the number of markers.

        Input contraints:

        - Must be sorted in ascending order.
        - Must be grouped based on chromosome. The chromosome boundary start
          and stop indices must be provided in ``chrgrp_stix`` and
          ``chrgrp_spix``, respectively.
    chrgrp_stix : numpy.ndarray
        Chromosome boundary start indices array of shape ``(c,)``.

        Where:

        - ``c`` is the number of chromosomes.
    chrgrp_spix : numpy.ndarray
        Chromosome boundary stop indices array of shape ``(c,)``.

        Where:

        - ``c`` is the number of chromosomes.

    Returns
    -------
    out : numpy.ndarray
        Array of shape ``(c,)`` containing the number of haplotype blocks
        assigned to each chromosome.

        Where:

        - ``c`` is the number of chromosomes.
    """
    # get number of chromosomes
    nchr = len(chrgrp_stix)

    # make sure that number of haplotype blocks >= number of chromosomes
    if nhaploblk < nchr:
        raise ValueError(
            "Number of haplotype blocks (nhaploblk = {0}) ".format(nhaploblk) +
            "is less than the number of chromosomes (nchr = {1})".format(nchr)
        )

    # calculate genetic lengths of each chromosome in Morgans/cM/other
    genlen = genpos[chrgrp_spix-1] - genpos[chrgrp_stix]

    # calculate ideal number of markers per chromosome
    nhaploblk_ideal = (nhaploblk / genlen.sum()) * genlen

    # calculate number of chromosome markers, assuming at least one per chromosome
    nhaploblk_chrom = numpy.ones(nchr, dtype = "int")    # start with min of one

    # greedily increment number of haplotype blocks assigned to chromosomes
    for i in range(nhaploblk - nchr):               # forces conformance to self.nhaploblk
        diff = nhaploblk_chrom - nhaploblk_ideal    # take actual - ideal
        ix = diff.argmin()                          # get index of lowest difference
        nhaploblk_chrom[ix] += 1                    # increment at lowest index

    return nhaploblk_chrom

def haplobin(
        nhaploblk_chrom: numpy.ndarray, 
        genpos: numpy.ndarray, 
        chrgrp_stix: numpy.ndarray, 
        chrgrp_spix: numpy.ndarray
    ) -> numpy.ndarray:
    """
    Given the number of haplotype blocks to give to each chromosome across the
    genome, assign bins for each marker

    Parameters
    ----------
    nhaploblk_chrom : numpy.ndarray
        Array of shape ``(c,)`` containing the total number of haplotype blocks
        to assign to each chromosome across the genome.

        Where:

        - ``c`` is the number of chromosomes.
    genpos : numpy.ndarray
        Array of shape ``(p,)`` containing genetic positions of markers across
        the genome.

        Where:

        - ``p`` is the number of markers.

        Input contraints:

        - Must be sorted in ascending order.
        - Must be grouped based on chromosome. The chromosome boundary start
          and stop indices must be provided in ``chrgrp_stix`` and
          ``chrgrp_spix``, respectively.
    chrgrp_stix : numpy.ndarray
        Chromosome boundary start indices array of shape ``(c,)``.

        Where:

        - ``c`` is the number of chromosomes.
    chrgrp_spix : numpy.ndarray
        Chromosome boundary stop indices array of shape ``(c,)``.

        Where:

        - ``c`` is the number of chromosomes.

    Returns
    -------
    out : numpy.ndarray
        Array of shape ``(p,)`` containing haplotype bin assignments for each
        marker across the genome. Bins are assigned starting at ``0`` and
        strictly increase as the index in the array increases due to constraints
        placed on ``genpos``.

        Where:

        - ``p`` is the number of markers.

        Example output::

            [0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3]
    """
    # initialize output array
    haplobin = numpy.empty(     # create empty array
        len(genpos),            # same length as number of markers
        dtype = "int"           # integer array
    )

    # get number of chromosomes
    nchr = len(chrgrp_stix)

    k = 0                           # haplotype block counter
    for i in range(nchr):           # for each chromosome
        nhap = nhaploblk_chrom[i]   # get number of haplotype blocks on chromosome
        stix = chrgrp_stix[i]       # get chromosome start index
        spix = chrgrp_spix[i]       # get chromosome stop index
        hbound = numpy.linspace(    # calculate haplotype boundaries
            genpos[stix],           # start genetic position on chromosome
            genpos[spix-1],         # stop genetic position on chromosome
            nhap+1                  # nhap + 1 for chromosome tips
        )
        for j in range(nhap):               # for each haplotype block
            chrmap = genpos[stix:spix]      # get chromosome genetic positions
            lmask = chrmap >= hbound[j]     # select all >= lower haplotype block bound
            umask = chrmap <= hbound[j+1]   # select all <= upper haplotype block bound
            mask = lmask & umask            # select all: lower <= markers <= upper
            haplobin[stix:spix][mask] = k   # set haplotype bin to 'k'
            k += 1                          # increment bin count

    return haplobin

def haplobin_bounds(
        haplobin: numpy.ndarray
    ) -> tuple:
    """
    Calculate haplotype bin boundaries and lengths.

    Parameters
    ----------
    haplobin : numpy.ndarray
        Array of shape ``(p,)`` containing haplotype bin assignments for each
        marker across the genome. Bins are assigned starting at ``0`` and
        strictly increase as the index in the array increases.

        Where:

        - ``p`` is the number of markers.

        Example input::

            [0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3]

    Returns
    -------
    out : tuple
        A tuple of length 3 containing ``(hstix, hspix, hlen)``.

        Where:

        - ``hstix`` is a ``numpy.ndarray`` containing haplotype bin boundary
          start indices.
        - ``hspix`` is a ``numpy.ndarray`` containing haplotype bin boundary
          stop indices.
        - ``hlen`` is a ``numpy.ndarray`` containing haplotype bin lengths.
    """
    hstix = [0]                         # starting indices
    hspix = []                          # stopping indices
    prev = haplobin[0]                  # get first group label
    for i in range(1, len(haplobin)):   # for each group label
        if haplobin[i] != prev:         # if the label is different
            hspix.append(i)             # add the stop index for prev
            prev = haplobin[i]          # get next label
            hstix.append(i)             # add the start index for next prev
    hspix.append(len(haplobin))         # append last stop index
    hstix = numpy.int_(hstix)           # convert to ndarray
    hspix = numpy.int_(hspix)           # convert to ndarray
    hlen = hspix - hstix                # get lengths of each haplotype group

    return hstix, hspix, hlen

def haplomat(
        nhaploblk: int,
        genomemat: numpy.ndarray,
        genpos: numpy.ndarray,
        chrgrp_stix: numpy.ndarray,
        chrgrp_spix: numpy.ndarray,
        chrgrp_len: numpy.ndarray,
        u_a: numpy.ndarray,
    ) -> numpy.ndarray:
    """
    Calculate a haplotype matrix from a genome matrix and model.

    Parameters
    ----------
    nhaploblk : Integral
        Number of haplotype blocks
    genomemat : numpy.ndarray
        A genome matrix of shape ``(m,n,p)``.

        Where:

        - ``m`` is the number of chromsome phases.
        - ``n`` is the number of individuals.
        - ``p`` is the number of markers.
    genpos : numpy.ndarray
        Array of shape ``(p,)`` containing genetic positions of markers across
        the genome.

        Where:

        - ``p`` is the number of markers.

        Input contraints:

        - Must be sorted in ascending order.
        - Must be grouped based on chromosome. The chromosome boundary start
          and stop indices must be provided in ``chrgrp_stix`` and
          ``chrgrp_spix``, respectively.
    chrgrp_stix : numpy.ndarray
        Chromosome boundary start indices array of shape ``(c,)``.

        Where:

        - ``c`` is the number of chromosomes.
    chrgrp_spix : numpy.ndarray
        Chromosome boundary stop indices array of shape ``(c,)``.

        Where:

        - ``c`` is the number of chromosomes.
    chrgrp_len : numpy.ndarray
        Chromosome length array of shape ``(c,)``.

        Where:

        - ``c`` is the number of chromosomes.
    u_a : numpy.ndarray
        Array of shape ``(p,t)`` containing additive marker effect estimates.

        Where:

        - ``p`` is the number of markers.
        - ``t`` is the number of traits.

    Returns
    -------
    out : numpy.ndarray
        A haplotype effect matrix of shape ``(m,n,b,t)``.

        Where:

        - ``m`` is the number of chromsome phases.
        - ``n`` is the number of individuals.
        - ``b`` is the number of haplotype blocks.
        - ``t`` is the number of traits.
    """
    if (chrgrp_stix is None) or (chrgrp_spix is None):
        raise RuntimeError("markers are not sorted by chromosome position")

    # get number of chromosomes
    nchr = len(chrgrp_stix)

    if nhaploblk < nchr:
        raise RuntimeError("number of haplotype blocks is less than the number of chromosomes")

    # calculate number of marker blocks to assign to each chromosome
    nblk = nhaploblk_chrom(nhaploblk, genpos, chrgrp_stix, chrgrp_spix)

    # ensure there are enough markers per chromosome
    if numpy.any(nblk > chrgrp_len):
        raise RuntimeError(
            "number of haplotype blocks assigned to a chromosome greater than number of available markers"
        )

    # calculate haplotype bins
    hbin = haplobin(nblk, genpos, chrgrp_stix, chrgrp_spix)

    # define shape
    s = (genomemat.shape[0], genomemat.shape[1], nhaploblk, u_a.shape[1]) # (m,n,b,t)

    # allocate haplotype matrix
    hmat = numpy.empty(s, dtype = u_a.dtype)   # (m,n,b,t)

    # get boundary indices
    hstix, hspix, hlen = haplobin_bounds(hbin)

    # OPTIMIZE: perhaps eliminate one loop using dot function
    # fill haplotype matrix
    for i in range(hmat.shape[3]):                          # for each trait
        for j,(st,sp) in enumerate(zip(hstix,hspix)):       # for each haplotype block
            hmat[:,:,j,i] = genomemat[:,:,st:sp].dot(u_a[st:sp,i])  # take dot product and fill

    return hmat
