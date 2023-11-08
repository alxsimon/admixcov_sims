#!python
import demes
import demesdraw
import msprime
# import math
import numpy as np
import stdpopsim
# import pickle

# import funcs as fn
# inputs
demes_file = snakemake.input['demes_file']
# outputs
trees_file = snakemake.output['trees_file']
# params
rec = snakemake.params['recombination_rate']
census_time = snakemake.params['census_time']
n_sample = snakemake.params['n_sample']
sampling_times = snakemake.params['sampling_times']

graph = demes.load(demes_file)

demography = msprime.Demography.from_demes(graph)
demography.add_census(time=census_time)
demography.sort_events()

# sampling
samples = []
for d in graph.demes:
	samples += [
		msprime.SampleSet(n_sample, population=d.name, time=t)
		for t in sampling_times
		if (t < d.epochs[0].start_time) & (t >= d.epochs[-1].end_time)
	]
	if (census_time < d.epochs[0].start_time) & (census_time >= d.epochs[-1].end_time):
		samples.append(msprime.SampleSet(n_sample, population=d.name, time=census_time))

# Contig setup
species = stdpopsim.get_species("HomSap")
# contigs = [
# 	species.get_contig(chr)
# 	for chr in ['chr22']
# ]

# Simulation
ts = msprime.sim_ancestry(
	samples=samples,
	ploidy=2,
	# recombination_rate=rate_map,
	sequence_length=1e8,
	recombination_rate=rec,
	demography=demography,
)

ts = msprime.sim_mutations(
	ts,
	rate=1e-8,
)

# drop sites with recurrent mutations
ts = ts.delete_sites(np.where([len(s.mutations) > 1 for s in ts.sites()])[0])

ts.dump(trees_file)
