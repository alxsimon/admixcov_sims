rule sim_msprime_simple_scenarios:
	input:
		demes_file = 'results/simulations/scenario_{sc}.yaml',
	output:
		trees_file = 'results/simulations/sim_msprime_scenario_{sc}.trees',
		model_plot = 'results/simulations/sim_scenario_{sc}.svg',
	params:
		census_time = 200,
		n_sample = 100,
		sampling_times = [200, 120, 100, 80, 60, 40, 20, 0],
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/sim_msprime_simple_scenarios.py'


rule sim_slim_sel_simple_scenarios:
	input:
		demes_file = 'results/simulations/scenario_{sc}.json',
	output:
		trees_file = 'results/simulations/sim_slim_sel_scenario_{sc}.trees',
		pheno_file = 'results/simulations/sim_slim_sel_scenario_{sc}_pheno.tsv',
	params:
		census_time = 200,
		n_sample = 100,
		sampling_times = 'c(200, 120, 100, 80, 60, 40, 20, 0)',
		shift_size = 1.0,
		shift_delay = 120, # delay of shift from admix_start
	log: 
		"logs/sim_slim_sel_simple_scenarios_{sc}.log"
	conda:
		"../envs/popgensim.yaml"
	shell:
		'''
		slim \
		-d 'JSON_FILE="{input.demes_file}"' \
		-d 'TREES_FILE="{output.trees_file}"' \
		-d 'PHENO_FILE="{output.pheno_file}"' \
		-d 'backward_sampling={params.sampling_times}' \
		-d 'N_sample={params.n_sample}' \
		-d 'census_time={params.census_time}' \
		-d 'shift_size={params.shift_size}' \
		-d 'shift_delay={params.shift_delay}' \
		workflow/scripts/sim_slim_sel_simple_scenarios.slim \
		> {log}
		'''


rule sim_msprime_simple_3pop_Patterson_like:
	input:
		demes_file = 'resources/model_simple_3pop_Patterson-like.yaml',
	output:
		trees_file = 'results/simulations/sim_3pop_simple_Patterson-like.trees',
		model_plot = 'results/simulations/sim_3pop_Patterson-like.svg',
	params:
		census_time = 200,
		n_sample = 100,
		sampling_times = [200, 135, 115, 95, 75, 55, 0],
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/sim_msprime_simple_scenarios.py'


rule sim_msprime_eu:
	input:
		demes_file = lambda w: variant_models[w.variant],
	output:
		trees_file = 'results/simulations/sim_eu_{variant}.trees',
		model_plot = 'results/simulations/sim_eu_{variant}.svg',
	params:
		census_time = 210,
		n_sample = 50,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/sim_msprime_eu.py'

# rule sim_slim_simple_scenarios_sel:
# 	input:
# 		demes_file = 'resources/scenario_{sc}.yaml',
# 	output:
# 		trees_file = 'results/simulations/sim_scenario_{sc}.trees',
# 		model_plot = 'results/simulations/sim_scenario_{sc}.svg',
# 		# rate_map_pickle = 'results/simulations/sim_scenario_{sc}_rate_map.pickle',
# 	params:
# 		census_time = 210,
# 		n_sample = 100,
# 		sampling_times = [210, 205, 185, 165, 145, 125, 105, 85, 65, 0],
# 	conda:
# 		"../envs/popgensim.yaml"
# 	shell:
		# """
		# ../scripts/sim_msprime_simple_scenarios.py
		# """

# rule sim_msprime_2pop:
# 	input:
# 		demes_file = 'resources/model_2pop_{ver}.yaml',
# 	output:
# 		trees_file = 'results/simulations/sim_2pop_{ver}.trees',
# 		model_plot = 'results/simulations/sim_2pop_{ver}.svg',
# 		rate_map_pickle = 'results/simulations/sim_2pop_{ver}_rate_map.pickle',
# 	params:
# 		census_time = 210,
# 		n_sample = 1000,
# 		sampling_times = [205, 185, 165, 145, 105, 85, 45, 0],
# 	conda:
# 		"../envs/popgensim.yaml"
# 	script:
# 		'../scripts/sim_msprime_simple.py'


# rule sim_msprime_simple_2pop:
# 	input:
# 		demes_file = 'resources/model_simple_2pop.yaml',
# 	output:
# 		trees_file = 'results/simulations/sim_2pop_simple.trees',
# 		model_plot = 'results/simulations/sim_2pop.svg',
# 		rate_map_pickle = 'results/simulations/sim_2pop_rate_map.pickle',
# 	params:
# 		census_time = 210,
# 		n_sample = 1000,
# 		# mutation_start_time = 205,
# 	conda:
# 		"../envs/popgensim.yaml"
# 	script:
# 		'../scripts/sim_msprime_simple.py'

# rule sim_msprime_simple_2pop_const_pulse:
# 	input:
# 		demes_file = 'resources/model_simple_2pop_const_pulse.yaml',
# 	output:
# 		trees_file = 'results/simulations/sim_2pop_const_pulse.trees',
# 		model_plot = 'results/simulations/sim_2pop_const_pulse.svg',
# 		rate_map_pickle = 'results/simulations/sim_2pop_const_pulse_rate_map.pickle',
# 	params:
# 		census_time = 210,
# 		n_sample = 1000,
# 		# mutation_start_time = 205,
# 	conda:
# 		"../envs/popgensim.yaml"
# 	script:
# 		'../scripts/sim_msprime_simple.py'

# rule sim_msprime_simple_2pop_no_pulse:
# 	input:
# 		demes_file = 'resources/model_simple_2pop_no_pulse.yaml',
# 	output:
# 		trees_file = 'results/simulations/sim_2pop_no_pulse.trees',
# 		model_plot = 'results/simulations/sim_2pop_no_pulse.svg',
# 		rate_map_pickle = 'results/simulations/sim_2pop_no_pulse_rate_map.pickle',
# 	params:
# 		census_time = 210,
# 		n_sample = 1000,
# 		# mutation_start_time = 205,
# 	conda:
# 		"../envs/popgensim.yaml"
# 	script:
# 		'../scripts/sim_msprime_simple.py'


# rule sim_msprime_simple_3pop:
# 	input:
# 		demes_file = 'resources/model_simple_3pop.yaml',
# 	output:
# 		trees_file = 'results/simulations/sim_3pop_simple.trees',
# 		model_plot = 'results/simulations/sim_3pop.svg',
# 		rate_map_pickle = 'results/simulations/sim_3pop_rate_map.pickle',
# 	params:
# 		census_time = 210,
# 		n_sample = 1000,
# 		# mutation_start_time = 205,
# 	conda:
# 		"../envs/popgensim.yaml"
# 	script:
# 		'../scripts/sim_msprime_simple.py'


# rule sim_msprime_simple_3pop_Patterson_like:
# 	input:
# 		demes_file = 'resources/model_simple_3pop_Patterson-like.yaml',
# 	output:
# 		trees_file = 'results/simulations/sim_3pop_simple_Patterson-like.trees',
# 		model_plot = 'results/simulations/sim_3pop_Patterson-like.svg',
# 	params:
# 		census_time = 200,
# 		n_sample = 100,
# 	conda:
# 		"../envs/popgensim.yaml"
# 	script:
# 		'../scripts/sim_msprime_3pop_Patterson-like.py'