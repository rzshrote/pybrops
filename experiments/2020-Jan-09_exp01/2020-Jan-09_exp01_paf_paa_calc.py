# the purpose of this experiment is to examine the correlation between allele
# availability and allele frequency
# Hypothesis: weakly correlated

import numpy, pandas, time, os, sys
# hack to append into our path the parent directory for this file
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.realpath(__file__)
    )))
)
# import our libraries
from pybropt import objfn

################################################################################
# define several variables
pldd_filename = "paf_topology_struct_unequalafreq.tsv"
seed = 192020
n_indiv = 100
n_sel = 10
n_loci = 10000
n_pts = 1000
cycles = 1
indices = numpy.array([i for i in range(n_indiv)])
t = 1.0 # target allele freq

# make variable to store points dataframe
pts_df = None

# if the data exists, load it, else generate it
if not os.path.exists(pldd_filename):
    # seed rng
    numpy.random.seed(seed)
    # generate allele proportions (induce structure in the population)
    # structure = numpy.repeat(0.1, n_indiv*n_loci).reshape(n_indiv,n_loci)
    structure = numpy.random.uniform(0.0,0.25,(n_indiv,n_loci))
    # structure = numpy.tile(
    #     numpy.random.uniform(0.0,1.0,n_loci),
    #     (n_indiv,1)
    # )
    # numpy.repeat(
    #     0.1,
    #     n_indiv*n_loci
    #     #(numpy.random.uniform(0.05,0.75,n_indiv//2),numpy.random.uniform(0.125,0.15,n_indiv//2)),
    #     #n_loci
    # ).reshape(n_indiv,n_loci)
    # allocate genotype matrix
    geno = numpy.empty((2,n_indiv,n_loci), dtype='uint8')
    # generate binary marker data
    geno[0,:,:] = numpy.random.binomial(1, structure, (1,n_indiv,n_loci)).astype('uint8')
    geno[1,:,:] = numpy.random.binomial(1, structure, (1,n_indiv,n_loci)).astype('uint8')
    # generate marker coefficients
    coeff = numpy.random.normal(0, 1, n_loci)
    # calculate marker weight coefficients
    wcoeff = numpy.absolute(coeff)
    # calculate tfreq
    tfreq = numpy.where(coeff >= 0, t, 1.0-t)

    # make empty arrays for storing scores
    pts = numpy.empty((n_pts,4), dtype='float64')

    # computational loop
    for i in range(n_pts):
        # generate random sample
        rsel = numpy.random.choice(indices, n_sel, replace=False)
        # score everything
        pts[i,0] = t
        pts[i,1] = objfn.paa(rsel, geno, wcoeff, tfreq)
        pts[i,2] = objfn.paf(rsel, geno, wcoeff)
        pts[i,3] = objfn.pafd(rsel, geno, wcoeff, tfreq)

    # make a dataframe
    pts_df = pandas.DataFrame(pts, columns=["tfreq","paa","paf","pafd"])

    # save dataframe
    pts_df.to_csv(pldd_filename, index=False, sep='\t')
