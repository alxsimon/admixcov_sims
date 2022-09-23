
variant_models = {
	'A': 'resources/model_europe.yaml',
	'B': 'resources/model_europe_noEHGpulse.yaml',
}

rule sim_msprime:
	input:
		demes_file = lambda w: variant_models[w.variant],
	output:
		trees_file = 'results/simulations/europe_sim_{variant}.trees',
		model_plot = 'results/simulations/europe_sim_{variant}.svg',
		rate_map_pickle = 'results/simulations/europe_sim_{variant}_rate_map.pickle',
	params:
		census_time = 210,
		n_sample = 50,
		mutation_start_time = 205,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/sim_msprime.py'