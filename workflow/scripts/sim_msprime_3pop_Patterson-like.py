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
model_plot = snakemake.output['model_plot']
# rate_map_pickle = snakemake.output['rate_map_pickle']
# params
census_time = snakemake.params['census_time']
n_sample = snakemake.params['n_sample']
# mutation_start_time = snakemake.params['mutation_start_time']


graph = demes.load(demes_file)
ax = demesdraw.tubes(graph, log_time=True)
ax.figure.savefig(model_plot)

demography = msprime.Demography.from_demes(graph)
demography.add_census(time=census_time)
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
sampling_times = [205, 135, 115, 95, 75, 55, 0]
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
	for chr in ['chr22'] # , 'chr2', 'chr3']
]


# def merge_RateMaps(list_contigs):
# 	r_break = math.log(2)
# 	list_RateMaps = [c.recombination_map.map for c in contigs]
# 	merged_pos = [list_RateMaps[0].position]
# 	pos_offset = merged_pos[0][-1] + 1
# 	merged_rate = [list_RateMaps[0].rate]
# 	for R in list_RateMaps[1:]:
# 		merged_pos.append(R.position + pos_offset)
# 		merged_rate.append(np.concatenate([[r_break], R.rate]))
# 		pos_offset = merged_pos[-1][-1] + 1
# 	merged_pos = np.concatenate(merged_pos)
# 	merged_rate = np.concatenate(merged_rate)
# 	return msprime.RateMap(position=merged_pos, rate=merged_rate)


# rate_map = merge_RateMaps(contigs)
# # save this rate map
# with open(rate_map_pickle, 'wb') as fw:
# 	pickle.dump(rate_map, fw)


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
