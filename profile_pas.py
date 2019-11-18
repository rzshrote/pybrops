import numpy
import cProfile
import timeit

# population architect selection
def pas(rslice, geno, coeff, tfreq, tfreq_edge, tld, tld_edge,
        cycles, d, lgroup_size):
    # calculate pairwise marker weights
    # make triangle matrix
    marker_score = numpy.triu(coeff * coeff[:,None])
    # scale to sum to 1
    marker_score /= marker_score.sum()

    # calculate probabilities for recombination
    # fill with 0.5 (for independent assortment)
    prec = numpy.empty((len(d), len(d)), dtype=numpy.float64)
    prec.fill(0.5)

    # stop indices
    dsp = numpy.cumsum(lgroup_size)

    # start indices
    dst = dsp - lgroup_size[0]

    # for each start, stop index
    for st,sp in zip(dst,dsp):
        # make mesh for that array region
        tmpmesh = numpy.meshgrid(d[st:sp], d[st:sp])
        # calculate Haldane distance for that region
        prec[st:sp,st:sp] = 0.5 * (1 - numpy.exp(
                        -2 * numpy.abs(tmpmesh[0] - tmpmesh[1])))

    # select rows, THIS MAKES A COPY OF THE ARRAY (PERFORMANCE ISSUE)
    phase = geno[:,rslice,:].reshape(len(rslice)*geno.shape[0], geno.shape[2])

    # get the dimensions of the matrix
    rows, columns = phase.shape

    # convert the number of sequences (rows) to float64
    nseq_f64 = numpy.float64(rows)

    # calculate allele sums (axis=0), divide by nseq_f64 to get probabilities
    pfreq = phase.sum(0) / nseq_f64

    # calculate D; multiply by decay over 'cycles' based on recombination rate
    D = ((numpy.matmul(phase.transpose(), phase) / nseq_f64) - \
          numpy.matmul(pfreq[:,None], pfreq[None,:])) * \
          (1 - prec)**cycles

    # calculate A = p(1-p)
    A = pfreq * (1 - pfreq)

    # calculate r^2 matrix
    r_sq = (D**2) / numpy.matmul(A, A.transpose())

    # determine allele availability
    allele_avail = numpy.where(
        tfreq >= 1.0,
        pfreq > 0.0,
        numpy.where(
            tfreq > 0.0,
            numpy.logical_and(pfreq > 0.0, pfreq < 1.0),
            pfreq < 1.0
        )
    )

    # calculate pairwise allele availability
    pair_avail = allele_avail * allele_avail[:,None]

    # calculate difference between target and population
    diff_tp = tfreq - pfreq

    # calculate difference between target and edge
    diff_te = tfreq - tfreq_edge

    # make a mesh for differences between target and population; this is a list
    mesh_tp = numpy.meshgrid(diff_tp, diff_tp)

    # make a mesh for differences between edge and population; this is a list
    mesh_te = numpy.meshgrid(diff_te, diff_te)

    # calculate pairwise distances between target and population freq. loci
    # calculate pairwise distances between the edge and population freq. loci
    # finally divide population distance from target by distance to edge
    # do 1 - dist_prime
    allele_score = 1 - (numpy.sqrt(mesh_tp[0]**2, mesh_tp[1]**2) /
                        numpy.sqrt(mesh_te[0]**2, mesh_te[1]**2))

    # calculate the difference between target LD and population LD
    # calculate the difference between target LD and edge LD
    # calculate LD score component
    # modify the ld_score based on ambiguous LD states (r_sq == 0)
    # TODO: replace this with an approx. equal due to floating point errors
    ld_score = numpy.where(
        r_sq == 0.0,
        allele_avail,
        1 - ((tld - r_sq) / (tld - tld_edge))
    )

    tmpdiag = (marker_score * allele_score * ld_score)

    print(tmpdiag)

    # marker_score is a triangle matrix, so it should zero out.
    # calculate product and return the sum of the triangle matrix.
    return tmpdiag.sum()



################################################################################
################################################################################
################################################################################


# generate binary marker data
n_lines = 20
n_markers = 1000
n_phases = 2
n_selected = 10

# seed random number
numpy.random.seed(111019)
# make markers
markers = numpy.random.binomial(1, 0.1, (n_phases,n_lines,n_markers))
# make effects
effects = numpy.random.normal(0,2,n_markers)

#rslice, geno, coeff, tfreq, tfreq_edge, tld, tld_edge,
#        cycles, d, lgroup_size

# calculate target frequencies: if effect >= 0, then tfreq = 1; else tfreq = 0
tfreq = numpy.where(effects >= 0.0, 1.0, 0.0)
# calculate the edge frequencies: if tfreq > 0.5, then edge = 0; else edge = 1
tfreq_edge = numpy.where(tfreq >= 0.5, 0.0, 1.0)
# set the target ld
tr_sq = 1.0
# calculate edge target ld
tr_sq_edge = 0.0 if tr_sq >= 0.5 else 1.0
# calculate abs(effects)
abs_effects = numpy.abs(effects)
# simulate a single chromosome of length 'n_markers'
gmap = numpy.sort(numpy.random.uniform(0.0,1.0,n_markers))
# make linkage group sizes
gmap_group_sizes = numpy.array([n_markers])
# set the number of breeding cycles
generations = 5
# select a random subset of individuals to score
selections = numpy.random.choice(numpy.arange(n_lines), n_selected, False)

exec1 = 'numpy.matmul(tfreq[:,None], tfreq[None,:])'
exec2 = 'tfreq[:,None] @ tfreq[None,:]'
exec1_list = list()
exec2_list = list()
for i in range(25):
    exec1_list.append(timeit.timeit(exec1, number=1000, globals=globals()))
    exec2_list.append(timeit.timeit(exec2, number=1000, globals=globals()))
with open("matrix_mult.txt", "w") as handle:
    for a,b in zip(exec1_list, exec2_list):
        handle.write("%s\t%s\n" % (a,b))

# profile our function
# cProfile.run('pas(selections, markers, abs_effects, tfreq, tfreq_edge, tr_sq,\
#             tr_sq_edge, generations, gmap, gmap_group_sizes)')
