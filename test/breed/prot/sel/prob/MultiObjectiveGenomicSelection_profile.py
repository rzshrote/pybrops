"""
Code to profile performance of the objective function for MOGS
"""

import numpy
import timeit

def method1(pfreq: numpy.ndarray, tfreq: numpy.ndarray) -> numpy.ndarray:
    # mask for whether the goal is to fix the minor allele
    # (p,t)
    tfreq_fix_minor = (tfreq <= 0.0)
    # mask for whether the goal is to fix the major allele
    # (p,t)
    tfreq_fix_major = (tfreq >= 1.0)
    # mask for whether the foal is to maintain heterozygosity
    # (p,t)
    tfreq_maintain_het = ~(tfreq_fix_minor | tfreq_fix_major)
    # mask for whether the major allele has been lost (minor allele is fixed)
    # (p,1)
    pfreq_major_islost = (pfreq <= 0.0)
    # mask for whether the minor allele has been lost (major allele is fixed)
    # (p,1)
    pfreq_minor_islost = (pfreq >= 1.0)
    # mask for whether heterozygosity has been lost (minor or major allele is fixed)
    # (p,1)
    pfreq_heter_islost = (pfreq_major_islost | pfreq_minor_islost)
    # calculate allele unavailability
    # (p,t)
    allele_unavail = numpy.where(
        tfreq_fix_major,        # (p,t) if target fixation of major allele
        pfreq_major_islost,     # (p,1) score True if major allele is lost
        numpy.where(
            tfreq_maintain_het, # (p,t) else if target heterozygosity maintenance
            pfreq_heter_islost, # (p,1) score True if heterozygosity is lost
            pfreq_minor_islost, # (p,1) else if target fixation of minor allele, score True if minor allele is lost
        )
    )
    return allele_unavail

def method2(pfreq: numpy.ndarray, tfreq: numpy.ndarray) -> numpy.ndarray:
    # mask for whether the goal is to fix the minor allele
    # (p,t)
    tfreq_fix_minor = (tfreq <= 0.0)
    # mask for whether the goal is to fix the major allele
    # (p,t)
    tfreq_fix_major = (tfreq >= 1.0)
    # mask for whether the foal is to maintain heterozygosity
    # (p,t)
    tfreq_fix_heter = ~(tfreq_fix_minor | tfreq_fix_major)
    # mask for whether the major allele has been lost (minor allele is fixed)
    # (p,1)
    pfreq_major_islost = (pfreq <= 0.0)
    # mask for whether the minor allele has been lost (major allele is fixed)
    # (p,1)
    pfreq_minor_islost = (pfreq >= 1.0)
    # mask for whether heterozygosity has been lost (minor or major allele is fixed)
    # (p,1)
    pfreq_heter_islost = (pfreq_major_islost | pfreq_minor_islost)
    # calculate minor, major, heterozygosity penalties
    # (p,t)
    minor_penalty = (tfreq_fix_minor & pfreq_minor_islost)
    major_penalty = (tfreq_fix_major & pfreq_major_islost)
    heter_penalty = (tfreq_fix_heter & pfreq_heter_islost)
    # calculate allele unavailability
    # (p,t)
    allele_unavail = minor_penalty | major_penalty | heter_penalty
    return allele_unavail

# construct pfreq and tfreq
pfreq_test = numpy.array([1.0,   1.0,   1.0,   0.5,   0.5,   0.5,   0.0,   0.0,   0.0  ])
tfreq_test = numpy.array([1.0,   0.5,   0.0,   1.0,   0.5,   0.0,   1.0,   0.5,   0.0  ])
truth_test = numpy.array([False, True,  True,  False, False, False, True,  True,  False]) # whether tfreq is unattainable

method1_test = method1(pfreq_test, tfreq_test)
method2_test = method2(pfreq_test, tfreq_test)

assert numpy.all(method1_test == truth_test)
assert numpy.all(method2_test == truth_test)

# create random vectors for timing purposes
nmarker = 1000
ntrait = 3

zeros = numpy.repeat(0.0, 300)
ones  = numpy.repeat(1.0, 300)
rands = numpy.random.random(400)

pfreq_time = numpy.concatenate([zeros, ones, rands])
tfreq_time = numpy.repeat(pfreq_time, 3)

numpy.random.shuffle(pfreq_time)
numpy.random.shuffle(tfreq_time)

pfreq_time = pfreq_time.reshape(nmarker, 1)
tfreq_time = tfreq_time.reshape(nmarker, ntrait)

method1_time = timeit.timeit("method1(pfreq_time, tfreq_time)", number = 500000, globals=globals())
# method 2 is superior
method2_time = timeit.timeit("method2(pfreq_time, tfreq_time)", number = 500000, globals=globals())

################################################################################
################################################################################
################################################################################

class MOGS1:
    def __init__(
            self,
            geno: numpy.ndarray,
            tfreq: numpy.ndarray,
            mkrwt: numpy.ndarray,
            ploidy: int,
            **kwargs: dict
        ) -> None:
        self.geno = geno
        self.tfreq = tfreq
        self.mkrwt = mkrwt
        self.ploidy = ploidy
    @property
    def geno(self) -> numpy.ndarray:
        return self._geno
    @geno.setter
    def geno(self, value: numpy.ndarray) -> None:
        self._geno = value
    @property
    def ploidy(self) -> int:
        return self._ploidy
    @ploidy.setter
    def ploidy(self, value: int) -> None:
        self._ploidy = value
    @property
    def mkrwt(self) -> numpy.ndarray:
        return self._mkrwt
    @mkrwt.setter
    def mkrwt(self, value: numpy.ndarray) -> None:
        self._mkrwt = value
    @property
    def tfreq(self) -> numpy.ndarray:
        return self._tfreq
    @tfreq.setter
    def tfreq(self, value: numpy.ndarray) -> None:
        self._tfreq = value
        self._tfreq_fix_minor = (self._tfreq <= 0.0)
        self._tfreq_fix_major = (self._tfreq >= 1.0)
        self._tfreq_fix_heter = ~(self._tfreq_fix_minor | self._tfreq_fix_major)
    @property
    def tfreq_fix_minor(self) -> numpy.ndarray:
        return self._tfreq_fix_minor
    @property
    def tfreq_fix_heter(self) -> numpy.ndarray:
        return self._tfreq_fix_heter
    @property
    def tfreq_fix_major(self) -> numpy.ndarray:
        return self._tfreq_fix_major
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        # calculate the allele frequency of the selected subset
        # (n,p)[(k,),:,None] -> (p,1)
        pfreq = (1.0 / (self.ploidy * len(x))) * self.geno[x,:,None].sum(0)
        # mask for whether the major allele has been lost (minor allele is fixed)
        # (p,1)
        pfreq_major_islost = (pfreq <= 0.0)
        # mask for whether the minor allele has been lost (major allele is fixed)
        # (p,1)
        pfreq_minor_islost = (pfreq >= 1.0)
        # mask for whether heterozygosity has been lost (minor or major allele is fixed)
        # (p,1)
        pfreq_heter_islost = (pfreq_major_islost | pfreq_minor_islost)
        # calculate minor, major, heterozygosity penalties
        # (p,t)
        minor_penalty = (self.tfreq_fix_minor & pfreq_minor_islost)
        major_penalty = (self.tfreq_fix_major & pfreq_major_islost)
        heter_penalty = (self.tfreq_fix_heter & pfreq_heter_islost)
        # calculate allele unavailability
        # (p,t)
        allele_unavail = minor_penalty | major_penalty | heter_penalty
        # calculate the allele unavailability
        # (p,t) * (p,t) -> (p,t)
        # (p,t) -> (t,)
        pau = (self.mkrwt * allele_unavail).sum(0)
        # calculate the manhattan distance and PAFD
        # (p,t) -> (t,)
        pafd = (self.mkrwt * numpy.absolute(self.tfreq - pfreq)).sum(0)
        # concatenate to make MOGS vector
        # (t,) and (t,) -> (t + t,)
        out = numpy.concatenate([pau, pafd])
        return out

class MOGS2:
    def __init__(
            self,
            geno: numpy.ndarray,
            tfreq: numpy.ndarray,
            mkrwt: numpy.ndarray,
            ploidy: int,
            **kwargs: dict
        ) -> None:
        self.geno = geno
        self.tfreq = tfreq
        self.mkrwt = mkrwt
        self.ploidy = ploidy
        self.cache = None
    @property
    def geno(self) -> numpy.ndarray:
        return self._geno
    @geno.setter
    def geno(self, value: numpy.ndarray) -> None:
        self._geno = value
    @property
    def ploidy(self) -> int:
        return self._ploidy
    @ploidy.setter
    def ploidy(self, value: int) -> None:
        self._ploidy = value
    @property
    def mkrwt(self) -> numpy.ndarray:
        return self._mkrwt
    @mkrwt.setter
    def mkrwt(self, value: numpy.ndarray) -> None:
        self._mkrwt = value
    @property
    def tfreq(self) -> numpy.ndarray:
        return self._tfreq
    @tfreq.setter
    def tfreq(self, value: numpy.ndarray) -> None:
        self._tfreq = value
        self._tfreq_fix_minor = (self._tfreq <= 0.0)
        self._tfreq_fix_major = (self._tfreq >= 1.0)
        self._tfreq_fix_heter = ~(self._tfreq_fix_minor | self._tfreq_fix_major)
    @property
    def tfreq_fix_minor(self) -> numpy.ndarray:
        return self._tfreq_fix_minor
    @property
    def tfreq_fix_heter(self) -> numpy.ndarray:
        return self._tfreq_fix_heter
    @property
    def tfreq_fix_major(self) -> numpy.ndarray:
        return self._tfreq_fix_major
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        # get pointer to genotypes
        geno = self.geno
        # get the number of individuals
        nselect = len(x)
        # get the number of markers
        nmarker = geno.shape[1]
        # if the cache has not been created or is the wrong size, allocate memory
        if self.cache is None or self.cache.shape != (nselect,nmarker):
            self.cache = numpy.empty((nselect,nmarker), dtype=geno.dtype)
        # get pointer to cache
        # k = len(x); p = nmarker
        # (k,p)
        cache = self.cache
        # subset rows into the cache, hopefully without spending time allocating memory
        geno.take(x, axis = 0, out = cache)
        # calculate the allele frequency of the selected subset
        # (n,p)[(k,),:,None] -> (p,1)
        pfreq = (1.0 / (self.ploidy * len(x))) * cache.sum(0, keepdims=True).T
        # mask for whether the major allele has been lost (minor allele is fixed)
        # (p,1)
        pfreq_major_islost = (pfreq <= 0.0)
        # mask for whether the minor allele has been lost (major allele is fixed)
        # (p,1)
        pfreq_minor_islost = (pfreq >= 1.0)
        # mask for whether heterozygosity has been lost (minor or major allele is fixed)
        # (p,1)
        pfreq_heter_islost = (pfreq_major_islost | pfreq_minor_islost)
        # calculate minor, major, heterozygosity penalties
        # (p,t)
        minor_penalty = (self.tfreq_fix_minor & pfreq_minor_islost)
        major_penalty = (self.tfreq_fix_major & pfreq_major_islost)
        heter_penalty = (self.tfreq_fix_heter & pfreq_heter_islost)
        # calculate allele unavailability
        # (p,t)
        allele_unavail = minor_penalty | major_penalty | heter_penalty
        # calculate the allele unavailability
        # (p,t) * (p,t) -> (p,t)
        # (p,t) -> (t,)
        pau = (self.mkrwt * allele_unavail).sum(0)
        # calculate the manhattan distance and PAFD
        # (p,t) -> (t,)
        pafd = (self.mkrwt * numpy.absolute(self.tfreq - pfreq)).sum(0)
        # concatenate to make MOGS vector
        # (t,) and (t,) -> (t + t,)
        out = numpy.concatenate([pau, pafd])
        return out

class MOGS3:
    def __init__(
            self,
            geno: numpy.ndarray,
            tfreq: numpy.ndarray,
            mkrwt: numpy.ndarray,
            ploidy: int,
            **kwargs: dict
        ) -> None:
        self.geno = geno
        self.tfreq = tfreq
        self.mkrwt = mkrwt
        self.ploidy = ploidy
        self.cache = None
    @property
    def geno(self) -> numpy.ndarray:
        return self._geno
    @geno.setter
    def geno(self, value: numpy.ndarray) -> None:
        self._geno = value
    @property
    def ploidy(self) -> int:
        return self._ploidy
    @ploidy.setter
    def ploidy(self, value: int) -> None:
        self._ploidy = value
    @property
    def mkrwt(self) -> numpy.ndarray:
        return self._mkrwt
    @mkrwt.setter
    def mkrwt(self, value: numpy.ndarray) -> None:
        self._mkrwt = value
    @property
    def tfreq(self) -> numpy.ndarray:
        return self._tfreq
    @tfreq.setter
    def tfreq(self, value: numpy.ndarray) -> None:
        self._tfreq = value
        self._tfreq_fix_minor = (self._tfreq <= 0.0)
        self._tfreq_fix_major = (self._tfreq >= 1.0)
        self._tfreq_fix_heter = ~(self._tfreq_fix_minor | self._tfreq_fix_major)
    @property
    def tfreq_fix_minor(self) -> numpy.ndarray:
        return self._tfreq_fix_minor
    @property
    def tfreq_fix_heter(self) -> numpy.ndarray:
        return self._tfreq_fix_heter
    @property
    def tfreq_fix_major(self) -> numpy.ndarray:
        return self._tfreq_fix_major
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        # get pointer to genotypes
        geno = self.geno
        # get the number of individuals
        nselect = len(x)
        # get the number of markers
        nmarker = geno.shape[1]
        # if the cache has not been created or is the wrong size, allocate memory
        if self.cache is None or self.cache.shape != (nselect,nmarker):
            self.cache = numpy.empty((nselect,nmarker), dtype=geno.dtype)
        # get pointer to cache
        # k = len(x); p = nmarker
        # (k,p)
        cache = self.cache
        # subset rows into the cache, hopefully without spending time allocating memory
        cache[:,:] = geno[x,:]
        # numpy.take(geno, x, axis = 0, out = cache)
        # calculate the allele frequency of the selected subset
        # (n,p)[(k,),:,None] -> (p,1)
        pfreq = (1.0 / (self.ploidy * len(x))) * cache.sum(0, keepdims=True).T
        # mask for whether the major allele has been lost (minor allele is fixed)
        # (p,1)
        pfreq_major_islost = (pfreq <= 0.0)
        # mask for whether the minor allele has been lost (major allele is fixed)
        # (p,1)
        pfreq_minor_islost = (pfreq >= 1.0)
        # mask for whether heterozygosity has been lost (minor or major allele is fixed)
        # (p,1)
        pfreq_heter_islost = (pfreq_major_islost | pfreq_minor_islost)
        # calculate minor, major, heterozygosity penalties
        # (p,t)
        minor_penalty = (self.tfreq_fix_minor & pfreq_minor_islost)
        major_penalty = (self.tfreq_fix_major & pfreq_major_islost)
        heter_penalty = (self.tfreq_fix_heter & pfreq_heter_islost)
        # calculate allele unavailability
        # (p,t)
        allele_unavail = minor_penalty | major_penalty | heter_penalty
        # calculate the allele unavailability
        # (p,t) * (p,t) -> (p,t)
        # (p,t) -> (t,)
        pau = (self.mkrwt * allele_unavail).sum(0)
        # calculate the manhattan distance and PAFD
        # (p,t) -> (t,)
        pafd = (self.mkrwt * numpy.absolute(self.tfreq - pfreq)).sum(0)
        # concatenate to make MOGS vector
        # (t,) and (t,) -> (t + t,)
        out = numpy.concatenate([pau, pafd])
        return out


# create random vectors for timing purposes
nmarker = 1000
ntrait = 3

zeros = numpy.repeat(0.0, 300)
ones  = numpy.repeat(1.0, 300)
rands = numpy.random.random(400)

tfreq = numpy.repeat(numpy.concatenate([zeros, ones, rands]), 3)
numpy.random.shuffle(tfreq)
tfreq = tfreq.reshape(nmarker, ntrait)
mkrwt = numpy.random.random((nmarker,ntrait))
geno = numpy.random.binomial(2, 0.1, (240,1000))

mogs1 = MOGS1(geno, tfreq, mkrwt, 2)
mogs2 = MOGS2(geno, tfreq, mkrwt, 2)
mogs3 = MOGS3(geno, tfreq, mkrwt, 2)

# test several solutions
# mogs2.latentfn(numpy.random.choice(240,40))

# test time of several random vectors
# MOGS1 is superior
mogs1_time = timeit.timeit("mogs1.latentfn(numpy.random.choice(240,40))", number = 100000, globals = globals())
mogs2_time = timeit.timeit("mogs2.latentfn(numpy.random.choice(240,40))", number = 100000, globals = globals())
mogs3_time = timeit.timeit("mogs3.latentfn(numpy.random.choice(240,40))", number = 100000, globals = globals())

# list times
(mogs1_time,mogs2_time,mogs3_time)
