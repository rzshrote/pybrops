import numpy
import pandas

n_trials = 120

selections_list = []

def read_tsv_drop_cols(fname):
    df = pandas.read_csv(
        fname,
        sep = '\t' #,
        # usecols = lambda x: x not in ["seed", "score", "gebv"]
    )
    return df

for i in range(n_trials):
    # write population to file
    icpso_sel_fname = "data/" + "opv_icpso_sel" + str(i+1).zfill(4) + ".tsv"
    hc_sa_set_sel_fname = "data/" + "opv_hc_sa_set_sel" + str(i+1).zfill(4) + ".tsv"

    selections_list.append(read_tsv_drop_cols(hc_sa_set_sel_fname))
    selections_list.append(read_tsv_drop_cols(icpso_sel_fname))

    print("Read trial %d" % i)

selections_df = pandas.concat(selections_list)

selections_df.to_csv(
    "sel_merge.tsv",
    sep = '\t',
    index = False
)
