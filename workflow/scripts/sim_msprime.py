#!python
import demes
import msprime
import math
import numpy as np
import stdpopsim

try: 
	snakemake
except:
	print("Not running in Snakemake")
	print("Defining variables manually")
	# inputs
	demes_file = '../../resources/model_europe.yaml'
	# outputs
	trees_file = '../../results/simulations/europe_sim_test.trees'
	# params
	census_time = 210
	n_sample = 50
	mutation_start_time = 205
else:
	# inputs
	demes_file = snakemake.input['demes_file']
	# outputs
	trees_file = snakemake.output['trees_file']
	# params
	census_time = snakemake.params['census_time']
	n_sample = snakemake.params['n_sample']
	mutation_start_time = snakemake.params['mutation_start_time']


graph = demes.load(demes_file)

demography = msprime.Demography.from_demes(graph)
demography.add_census(time=census_time)
demography.sort_events()

# sampling
sampling_times = [205, 195, 185, 145, 135, 125, 100, 90, 80, 70, 0]
samples = [
	msprime.SampleSet(n_sample, population='Pop0', time=t)
	for t in sampling_times
]

# Contig setup
species = stdpopsim.get_species("HomSap")
contigs = [
	species.get_contig(
		chr,
		genetic_map='HapMapII_GRCh37',
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

# Simulation
ts = msprime.sim_ancestry(
	samples=samples,
	ploidy=2,
	recombination_rate=rate_map,
	demography=demography,
	# random_seed=,
)

ts = msprime.sim_mutations(
	ts,
	rate=contigs[0].mutation_rate,
	# keep=True, # keep existing mutations
	start_time=mutation_start_time,
	# random_seed=,
)

ts.dump(trees_file)