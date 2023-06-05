rule sim_msprime_simple_scenarios:
	input:
		demes_file = 'results_test/simulations/scenario_{sc}.yaml',
	output:
		trees_file = 'results_test/simulations/sim_msprime_rep/sim_msprime_scenario_{sc}_{rep}.trees',
	params:
		census_time = 200,
		n_sample = 100,
		sampling_times = [200, 160, 140, 120, 100, 80, 60, 40, 20, 0],
	resources:
		mem_mb = 5_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/sim_msprime_simple_scenarios.py'


rule sim_slim_sel_simple_scenarios:
	input:
		demes_file = 'results_test/simulations/scenario_{sc}.json',
	output:
		trees_file = temp('results_test/simulations/sim_slim_sel_rep/raw_sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}_{rep}.trees'),
		pheno_file = 'results_test/simulations/sim_slim_sel_rep/sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}_{rep}_pheno.tsv',
	params:
		census_time = 200,
		n_sample = 100,
		sampling_times = 'c(200, 160, 140, 120, 100, 80, 60, 40, 20, 0)',
		shift_delay = lambda w: 200 - int(w.time), # delay of shift from admix_start
	resources:
		mem_mb = 5_000,
	log: 
		"logs/sim_slim_sel_simple_scenarios_{sc}_{type}_t{time}_s{ssize}_{rep}.log"
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
		-d 'shift_type="{wildcards.type}"' \
		-d 'shift_size={wildcards.ssize}' \
		-d 'shift_delay={params.shift_delay}' \
		workflow/scripts/sim_slim_sel_simple_scenarios.slim \
		> {log}
		'''


rule sim_slim_sel_postprocessing:
	input:
		trees_file = 'results_test/simulations/sim_slim_sel_rep/raw_sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}_{rep}.trees',
		demes_file = 'results_test/simulations/scenario_{sc}.json',
	output:
		trees_file = 'results_test/simulations/sim_slim_sel_rep/sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}_{rep}.trees',
	params:
		neutral_mut_rate = 1.29e-08, # stdpopsim Human mutation rate
	resources:
		mem_mb = 5_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/sim_slim_postprocessing.py'


rule sim_msprime_europe_Patterson2022:
	input:
		demes_file = 'resources/model_europe_Patterson2022.yaml',
	output:
		trees_file = 'results_test/simulations/sim_msprime_europe_Patterson2022/sim_msprime_europe_Patterson2022_{rep}.trees',
	params:
		census_time = 200,
		n_sample = 300,
	resources:
		mem_mb = 5_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/sim_msprime_europe_Patterson2022.py'

