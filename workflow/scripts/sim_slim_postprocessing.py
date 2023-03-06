import tskit
import demes
import pyslim
import msprime
import numpy as np

in_trees = snakemake.input['trees_file']
demes_file = snakemake.input['demes_file']
out_trees = snakemake.output['trees_file']
mut_rate = snakemake.params['neutral_mut_rate']

ts = pyslim.update( # as we use SLiM v3.7
    tskit.load(in_trees)
)
graph = demes.load(demes_file)

Ne = graph.demes[0].epochs[0].start_size

# recap and add neutral mutations
ts_recap = pyslim.recapitate(ts, ancestral_Ne=Ne)
ts_mut = msprime.sim_mutations(
    ts_recap,
    rate=mut_rate,
    keep=False, # discard existing mutations
)
ts_mut = ts_mut.delete_sites(
    np.where([
        len(s.mutations) > 1 for s in ts_mut.sites()
    ])[0]
)

ts_mut.dump(out_trees)