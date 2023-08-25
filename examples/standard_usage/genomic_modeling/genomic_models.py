#!/usr/bin/env python3

import numpy

# import GenomicModel classes (abstract interface classes)
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.gmod.NonlinearGenomicModel import NonlinearGenomicModel
from pybrops.model.gmod.LinearGenomicModel import LinearGenomicModel
from pybrops.model.gmod.CoancestryLinearGenomicModel import CoancestryLinearGenomicModel
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.AdditiveDominanceLinearGenomicModel import AdditiveDominanceLinearGenomicModel
from pybrops.model.gmod.AdditiveDominanceEpistaticLinearGenomicModel import AdditiveDominanceEpistaticLinearGenomicModel

# import dense genomic models (concrete implementation classes)
from pybrops.model.gmod.DenseLinearGenomicModel import DenseLinearGenomicModel
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix

###
### Genomic model object construction
###

#
# Construction from NumPy arrays
#

# model parameters
nfixed = 1      # number of fixed effects
ntrait = 2      # number of traits
nmisc = 0       # number of miscellaneous random effects
nadditive = 50  # number of additive marker effects

# create dummy values
beta = numpy.random.random((nfixed,ntrait))
u_misc = numpy.random.random((nmisc,ntrait))
u_a = numpy.random.random((nadditive,ntrait))
trait = numpy.array(
    ["Trait"+str(i+1).zfill(2) for i in range(ntrait)],
    dtype = object
)

# create additive linear genomic model
algmod = DenseAdditiveLinearGenomicModel(
    beta = beta,
    u_misc = u_misc,
    u_a = u_a,
    trait = trait,
    model_name = "example",
    params = None
)

###
### Construct dummy genotypes
###

# shape parameters for random genotypes
ntaxa = 100
nvrnt = nadditive
ngroup = 20
nchrom = 10
ploidy = 2

# create random genotypes
mat = numpy.random.randint(0, ploidy+1, size = (ntaxa,nvrnt)).astype("int8")

# create taxa names
taxa = numpy.array(
    ["taxon"+str(i+1).zfill(3) for i in range(ntaxa)], 
    dtype = object
)

# create taxa groups
taxa_grp = numpy.random.randint(1, ngroup+1, ntaxa)
taxa_grp.sort()

# create marker variant chromsome assignments
vrnt_chrgrp = numpy.random.randint(1, nchrom+1, nvrnt)
vrnt_chrgrp.sort()

# create marker physical positions
vrnt_phypos = numpy.random.choice(1000000, size = nvrnt, replace = False)
vrnt_phypos.sort()

# create marker variant names
vrnt_name = numpy.array(
    ["SNP"+str(i+1).zfill(4) for i in range(nvrnt)],
    dtype = object
)

# create a genotype matrix from scratch using NumPy arrays
gmat = DenseGenotypeMatrix(
    mat = mat,
    taxa = taxa,
    taxa_grp = taxa_grp, 
    vrnt_chrgrp = vrnt_chrgrp,
    vrnt_phypos = vrnt_phypos, 
    vrnt_name = vrnt_name, 
    vrnt_genpos = None,
    vrnt_xoprob = None, 
    vrnt_hapgrp = None, 
    vrnt_hapalt = None,
    vrnt_hapref = None, 
    vrnt_mask = None,
    ploidy = ploidy
)


# general properties
algmod.model_name
algmod.params

# model properties
algmod.beta
algmod.u
algmod.u_a
algmod.u_misc

# trait properties
algmod.ntrait
algmod.trait

# methods

###
### Model fitting methods
###

algmod.fit_numpy()
algmod.fit()

algmod.predict_numpy()
algmod.predict()

algmod.score_numpy()
algmod.score()

algmod.gebv_numpy()
algmod.gebv()

algmod.var_G_numpy()
algmod.var_G()
algmod.var_A_numpy()
algmod.var_A()
algmod.var_a_numpy()
algmod.var_a()
algmod.bulmer_numpy()
algmod.bulmer()

algmod.usl_numpy()
algmod.usl()
algmod.lsl_numpy()
algmod.lsl()

algmod.facount()
algmod.fafreq()
algmod.faavail()
algmod.faavailval(gmat)
algmod.fafixed()
algmod.fafixedval()

algmod.dacount()
algmod.dafreq()
algmod.daavail()
algmod.daavailval()
algmod.dafixed()
algmod.dafixedval()

algmod.polyval()

algmod.to_hdf5()
algmod.from_hdf5()