import glob
import pandas
import numpy

pop_files = "*_pop*.tsv"
sel_files = "*_sel*.tsv"

pop_files_list = [
    pandas.read_csv(
        file,
        sep='\t',
        header=0
    ) for file in glob.glob(pop_files)
]

sel_files_list = [
    pandas.read_csv(
        file,
        sep='\t',
        header=0
    ) for file in glob.glob(sel_files)
]

pop_files_df = pandas.concat(
    pop_files_list,
    axis=0
)

sel_files_df = pandas.concat(
    sel_files_list,
    axis=0
)

pop_files_df.to_csv(
    "pop.tsv",
    sep='\t',
    index=False
)

sel_files_df.to_csv(
    "sel.tsv",
    sep='\t',
    index=False
)
