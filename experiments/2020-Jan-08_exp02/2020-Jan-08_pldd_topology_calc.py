# The purpose of this experiment is to examine the relationship between
# paa, pafd and pa as they vary with target frequency

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
pldd_filename = "pldd_topology.tsv"
seed = 182020
n_indiv = 100
n_sel = 10
n_loci = 500
n_pts = 1000 # don't want this too high due to LD matrix math complexity O(n^3)
fav_tfreq = numpy.repeat(1.0, n_pts) # numpy.linspace(0.0,1.0,n_pts)
tld = numpy.repeat(1.0, n_pts) # numpy.linspace(0.0,1.0,n_pts)
cycles = 1
indices = numpy.array([i for i in range(n_indiv)])


# make variable to store points dataframe
pts_df = None

# if the data exists, load it, else generate it
if not os.path.exists(pldd_filename):
    # seed rng
    numpy.random.seed(seed)
    # generate allele proportions (induce structure in the population)
    structure = numpy.random.uniform(0.0,0.2,(n_indiv,n_loci))
    # numpy.tile(
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
    # make target allele frequency variable
    tfreq = None
    # make genetic maps
    gmap = numpy.random.uniform(0.0, 1.5, n_loci)
    for b in range(0, n_loci, n_loci//10):
        gmap[b:b+(n_loci//10)] = numpy.sort(gmap[b:b+(n_loci//10)])
    lgroup = numpy.repeat(n_loci//10, 10)

    # make empty arrays for storing scores
    pts = numpy.empty((n_pts,6), dtype='float64')

    # computational loop
    for i,t,d in zip(range(n_pts),fav_tfreq,tld):
        # generate random sample
        rsel = numpy.random.choice(indices, n_sel, replace=False)
        # recalculate tfreq
        tfreq = numpy.where(coeff >= 0, t, 1.0-t)
        # score everything
        pts[i,0] = t
        pts[i,1] = d
        pts[i,2] = objfn.paa(rsel, geno, wcoeff, tfreq)
        pts[i,3] = objfn.pafd(rsel, geno, wcoeff, tfreq)
        pts[i,4] = objfn.pafd_prime(rsel, geno, wcoeff, tfreq)
        pts[i,5] = objfn.pldd(rsel, geno, wcoeff, tfreq, d, gmap, lgroup, cycles)

    # make a dataframe
    pts_df = pandas.DataFrame(pts, columns=["tfreq","tld","paa","pafd","pafd_prime","pldd"])

    # save dataframe
    pts_df.to_csv(pldd_filename, index=False, sep='\t')
