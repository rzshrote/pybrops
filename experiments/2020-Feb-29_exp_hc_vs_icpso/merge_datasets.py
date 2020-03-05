import numpy
import pandas

n_trials = 120

population_list = []
selections_list = []

def read_csv_add_cols(fname, trial):
    pop_df = pandas.read_csv(
        fname,
        sep = '\t',
        usecols = ["method", "algorithm", "seed", "cycle", "score", "gebv"]
    )
    # get mean gebvs for each cycle
    mean_gebv = []
    median_gebv = []
    min_gebv = []
    max_gebv = []
    std_gebv = []
    for name,group in pop_df.groupby("cycle"):
        tmp = group["gebv"].values
        mean_gebv.append(numpy.repeat(tmp.mean(), len(tmp)))
        median_gebv.append(numpy.repeat(numpy.median(tmp), len(tmp)))
        min_gebv.append(numpy.repeat(tmp.min(), len(tmp)))
        max_gebv.append(numpy.repeat(tmp.max(), len(tmp)))
        std_gebv.append(numpy.repeat(tmp.std(), len(tmp)))

    mean_gebv = numpy.concatenate(mean_gebv)
    median_gebv = numpy.concatenate(median_gebv)
    min_gebv = numpy.concatenate(min_gebv)
    max_gebv = numpy.concatenate(max_gebv)
    std_gebv = numpy.concatenate(std_gebv)

    pop_df["min_gebv"] = min_gebv
    pop_df["median_gebv"] = median_gebv
    pop_df["max_gebv"] = max_gebv
    pop_df["mean_gebv"] = mean_gebv
    pop_df["std_gebv"] = std_gebv
    pop_df["trial"] = trial

    return pop_df

for i in range(n_trials):
    # write population to file
    icpso_pop_fname = "data/" + "opv_icpso_pop" + str(i+1).zfill(4) + ".tsv"
    icpso_sel_fname = "data/" + "opv_icpso_sel" + str(i+1).zfill(4) + ".tsv"
    hc_sa_set_pop_fname = "data/" + "opv_hc_sa_set_pop" + str(i+1).zfill(4) + ".tsv"
    hc_sa_set_sel_fname = "data/" + "opv_hc_sa_set_sel" + str(i+1).zfill(4) + ".tsv"

    population_list.append(read_csv_add_cols(hc_sa_set_pop_fname, i+1))
    population_list.append(read_csv_add_cols(icpso_pop_fname, i+1))

    selections_list.append(read_csv_add_cols(hc_sa_set_sel_fname, i+1))
    selections_list.append(read_csv_add_cols(icpso_sel_fname, i+1))

    print("Read trial %d" % i)

population_df = pandas.concat(population_list)
selections_df = pandas.concat(selections_list)

population_df.to_csv(
    "pop_summary.tsv",
    sep = '\t',
    index = False
)
selections_df.to_csv(
    "sel_summary.tsv",
    sep = '\t',
    index = False
)
