#!python
import demes
# import demesdraw
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
# model_plot = snakemake.output['model_plot']
# rate_map_pickle = snakemake.output['rate_map_pickle']
# params
# census_time = snakemake.params['census_time']
n_sample = snakemake.params['n_sample']


graph = demes.load(demes_file)
# ax = demesdraw.tubes(graph, log_time=True)
# ax.figure.savefig(model_plot)

census_times = {'WHG': 200, 'ANA': 200, 'YAM': 150}

demography = msprime.Demography.from_demes(graph)
for ct in np.unique(list(census_times.values())):
	demography.add_census(time=ct)
demography.sort_events()

# cohorts = [
# 	'England.and.Wales_N',
# 	'England.and.Wales_C.EBA',
# 	'England.and.Wales_MBA',
# 	'England.and.Wales_LBA',
# 	'England.and.Wales_IA',
# 	'England.and.Wales_PostIA',
# 	'England.and.Wales_Modern',
# ]

# sampling
sampling_times = [150, 130, 110, 90, 70, 50, 0]
samples = []
for d in graph.demes:
	if d.name in ['NEO', 'WHG', 'ANA', 'YAM']:
		samples += [
			msprime.SampleSet(n_sample, population=d.name, time=t)
			for t in sampling_times
			if (t < d.epochs[0].start_time) & (t >= d.epochs[-1].end_time)
		]
for pop, ct in census_times.items():
	samples.append(msprime.SampleSet(n_sample, population=pop, time=ct))


# Contig setup
species = stdpopsim.get_species("HomSap")
contigs = [
	species.get_contig(chr)
	for chr in ['chr1']
]

# Simulation
ts = msprime.sim_ancestry(
	samples=samples,
	ploidy=2,
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
