import numpy
import pandas

################################################################################
# define several constants
seed = 112119
n_indiv = 200
n_loci = 10000
n_sel = 10
n_gen = 5
n_trials = 120

# generate binary marker data
numpy.random.seed(seed)
markers = numpy.empty((2,n_indiv,n_loci), dtype=numpy.uint8)
markers[0,:,:] = numpy.random.binomial(1, 0.1, (1,n_indiv,n_loci))
markers[1,:,:] = markers[0,:,:]
effects = numpy.random.normal(0, 1, n_loci)
pos_set = numpy.arange(n_indiv, dtype='int64')
pos_state = numpy.tile(pos_set, n_sel) # select 10
dim_sizes = numpy.repeat(n_indiv, n_sel)

# make genetic maps
gmap = numpy.random.uniform(0.0, 1.5, n_loci)
for b in range(0, n_loci, n_loci//10):
    gmap[b:b+(n_loci//10)] = numpy.sort(gmap[b:b+(n_loci//10)])
gmap_size = numpy.repeat(n_loci//10, 10)

marker_id_arr = numpy.arange(1, n_loci+1)
marker_effect_arr = effects.copy()
gmap_pos_arr = gmap.copy()
chr_arr = numpy.repeat(numpy.arange(1,11), n_loci//10)
################################################################################

df = pandas.read_csv(
    "sel_merge.tsv",
    sep = '\t',
    usecols = lambda x: x not in ["seed", "score", "gebv"]
)

groups = df.groupby(by=["algorithm", "cycle"])

algorithm = []
cycle = []
marker_id = []
marker_effect = []
chr = []
gmap_pos = []
marker_count = []

for name, group in groups:
    # gather markers
    tmp = numpy.fromiter(
        (group[m].values.sum() for m in group.columns[3:]),
        dtype = 'int64'
    )
    tmp = tmp[:len(tmp)//2] + tmp[len(tmp)//2:]

    algorithm.append(numpy.object_(numpy.repeat(name[0], len(tmp))))
    cycle.append(numpy.repeat(name[1], len(tmp)))
    marker_id.append(marker_id_arr)
    chr.append(chr_arr)
    gmap_pos.append(gmap_pos_arr)
    marker_effect.append(marker_effect_arr)
    marker_count.append(tmp)

# gather into dictionary
out_dict = {
    "algorithm" : numpy.concatenate(algorithm),
    "cycle" : numpy.concatenate(cycle),
    "marker_id" : numpy.concatenate(marker_id),
    "chr" : numpy.concatenate(chr),
    "gmap_pos" : numpy.concatenate(gmap_pos),
    "marker_effect" : numpy.concatenate(marker_effect),
    "marker_count" : numpy.concatenate(marker_count)
}

out_df = pandas.DataFrame(out_dict)

out_df.to_csv(
    "sel_allele_counts.tsv",
    sep = '\t',
    index = False
)
