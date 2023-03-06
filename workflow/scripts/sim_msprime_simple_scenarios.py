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
census_time = snakemake.params['census_time']
n_sample = snakemake.params['n_sample']
sampling_times = snakemake.params['sampling_times']
# mutation_start_time = snakemake.params['mutation_start_time']


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
contigs = [
	species.get_contig(chr)
	for chr in ['chr22'] # ['chr1'] # , 'chr2', 'chr3']
]

# Simulation
ts = msprime.sim_ancestry(
	samples=samples,
	ploidy=2,
	# recombination_rate=rate_map,
	recombination_rate=contigs[0].recombination_map,
	demography=demography,
)

ts = msprime.sim_mutations(
	ts,
	rate=contigs[0].mutation_rate,
)

# drop sites with recurrent mutations
ts = ts.delete_sites(np.where([len(s.mutations) > 1 for s in ts.sites()])[0])

ts.dump(trees_file)
