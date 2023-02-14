#!python
import demes
import demesdraw
import msprime
import math
import numpy as np
import stdpopsim
import pickle

try: 
	snakemake
except:
	print("Not running in Snakemake")
	print("Defining variables manually")
	import funcs as fn
	# inputs
	demes_file = 'resources/model_europe.yaml'
	# outputs
	trees_file = 'results/simulations/sim_eu_test.trees'
	model_plot = 'results/simulations/sim_eu_test.svg'
	rate_map_pickle = 'results/simulations/sim_eu_rate_map_test.pickle'
	# params
	census_time = 210
	n_sample = 50
	# mutation_start_time = 205
else:
	import funcs as fn
	# inputs
	demes_file = snakemake.input['demes_file']
	# outputs
	trees_file = snakemake.output['trees_file']
	model_plot = snakemake.output['model_plot']
	rate_map_pickle = snakemake.output['rate_map_pickle']
	# params
	census_time = snakemake.params['census_time']
	n_sample = snakemake.params['n_sample']
	# mutation_start_time = snakemake.params['mutation_start_time']


graph = demes.load(demes_file)
ax = demesdraw.tubes(graph, log_time=True)
ax.figure.savefig(model_plot)

demography = msprime.Demography.from_demes(graph)
demography.add_census(time=census_time)
demography.add_census(time=175)
demography.sort_events()

# sampling
sampling_times = [205, 195, 185, 175, 145, 135, 125, 100, 90, 80, 70, 0]
samples = [
	msprime.SampleSet(n_sample, population='Pop0', time=t)
	for t in sampling_times
]
samples.append(msprime.SampleSet(n_sample, population='Pop0', time=210))
samples.append(msprime.SampleSet(n_sample, population='Pop1', time=175))
samples.append(msprime.SampleSet(n_sample, population='Pop2', time=210))
samples.append(msprime.SampleSet(n_sample, population='Pop3', time=210))

# Contig setup
species = stdpopsim.get_species("HomSap")
# contig = species.get_contig('chr1', genetic_map='HapMapII_GRCh37')
# contig = species.get_contig('chr1', genetic_map='DeCodeSexAveraged_GRCh36')
contigs = [
	species.get_contig(
		chr,
		# length_multiplier=1,
		# genetic_map='HapMapII_GRCh37',
	)
	for chr in ['chr1', 'chr2', 'chr3']
]


def merge_RateMaps(list_contigs):
	r_break = math.log(2)
	list_RateMaps = [c.recombination_map.map for c in contigs]
	merged_pos = [list_RateMaps[0].position]
	pos_offset = merged_pos[0][-1] + 1
	merged_rate = [list_RateMaps[0].rate]
	for R in list_RateMaps[1:]:
		merged_pos.append(R.position + pos_offset)
		merged_rate.append(np.concatenate([[r_break], R.rate]))
		pos_offset = merged_pos[-1][-1] + 1
	merged_pos = np.concatenate(merged_pos)
	merged_rate = np.concatenate(merged_rate)
	return msprime.RateMap(position=merged_pos, rate=merged_rate)


rate_map = merge_RateMaps(contigs)
# rate_map = contig.recombination_map.map
# save this rate map
with open(rate_map_pickle, 'wb') as fw:
	pickle.dump(rate_map, fw)


# Simulation
ts = msprime.sim_ancestry(
	samples=samples,
	ploidy=2,
	recombination_rate=rate_map,
	# sequence_length=contig.recombination_map.get_length(),
	# recombination_rate=1.14856e-08,
	demography=demography,
	# random_seed=,
)

ts = msprime.sim_mutations(
	ts,
	rate=contigs[0].mutation_rate,
	# rate=contig.mutation_rate,
	# keep=True, # keep existing mutations
	# start_time=mutation_start_time,
	# random_seed=,
)

# drop sites with recurrent mutations
ts = ts.delete_sites(np.where([len(s.mutations) > 1 for s in ts.sites()])[0])

ts.dump(trees_file)
