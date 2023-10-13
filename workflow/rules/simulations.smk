rule sim_msprime_simple_scenarios:
	input:
		demes_file = 'results/simulations/scenario_{sc}.yaml',
	output:
		trees_file = 'results/simulations/sim_msprime_rep/sim_msprime_scenario_{sc}_{rep}.trees',
	params:
		census_time = 200,
		n_sample = 100,
		sampling_times = [200, 160, 140, 120, 100, 80, 60, 40, 20, 0],
	resources:
		mem_mb = 3_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/sim_msprime_simple_scenarios.py'


rule sim_slim_sel_simple_scenarios:
	input:
		demes_file = 'results/simulations/scenario_{sc}.json',
	output:
		trees_file = temp('results/simulations/sim_slim_sel_rep/raw_sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}_{rep}.trees'),
		pheno_file = 'results/simulations/sim_slim_sel_rep/sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}_{rep}_pheno.tsv',
	params:
		census_time = 200,
		n_sample = 100,
		sampling_times = 'c(200, 160, 140, 120, 100, 80, 60, 40, 20, 0)',
		shift_delay = lambda w: 200 - int(w.time), # delay of shift from admix_start
	resources:
		mem_mb = 9_000,
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
		trees_file = 'results/simulations/sim_slim_sel_rep/raw_sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}_{rep}.trees',
		demes_file = 'results/simulations/scenario_{sc}.json',
	output:
		trees_file = 'results/simulations/sim_slim_sel_rep/sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}_{rep}.trees',
	params:
		neutral_mut_rate = 1e-08,
	resources:
		mem_mb = 3_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/sim_slim_postprocessing.py'


rule sim_msprime_europe_uk:
	input:
		demes_file = 'resources/AncientEurope_4A21_mod.yaml',
	output:
		trees_file = 'results/simulations/sim_msprime_europe_uk/sim_msprime_europe_uk_{rep}.trees',
	params:
		n_sample = 300,
	resources:
		mem_mb = 9_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/sim_msprime_europe_uk.py'



rule sim_slim_sel_variable_high_sampling_freq:
	input:
		demes_file = 'results/simulations/scenario_{sc}.json',
	output:
		trees_file = 'results/simulations/sim_slim_sel_rep/raw_sim_slim_sel_inter_{sc}_{type}_t{time}_s{ssize}_r{rec}_{rep}.trees',
		pheno_file = 'results/simulations/sim_slim_sel_rep/sim_slim_sel_inter_{sc}_{type}_t{time}_s{ssize}_r{rec}_{rep}_pheno.tsv',
	params:
		census_time = 200,
		n_sample = 50,
		sampling_times = lambda w: f"c({','.join([str(x) for x in list(range(0, 201, 1))[::-1]])})",
		shift_delay = lambda w: 200 - int(w.time), # delay of shift from admix_start
	resources:
		mem_mb = 9_000,
	log:
		"logs/sim_slim_sel_high_interval_sampling_{sc}_{type}_t{time}_s{ssize}_r{rec}_{rep}.log"
	conda:
		"../envs/popgensim.yaml"
	shell:
		'''
		slim \
		-d 'JSON_FILE="{input.demes_file}"' \
		-d 'TREES_FILE="{output.trees_file}"' \
		-d 'PHENO_FILE="{output.pheno_file}"' \
		-d 'rec={wildcards.rec}'
		-d 'backward_sampling={params.sampling_times}' \
		-d 'N_sample={params.n_sample}' \
		-d 'census_time={params.census_time}' \
		-d 'shift_type="{wildcards.type}"' \
		-d 'shift_size={wildcards.ssize}' \
		-d 'shift_delay={params.shift_delay}' \
		workflow/scripts/sim_slim_sel_simple_scenarios.slim \
		> {log}
		'''

rule sim_slim_sel_variable_interval_postprocessing:
	input:
		trees_file = temp('results/simulations/sim_slim_sel_rep/raw_sim_slim_sel_inter_{sc}_{type}_t{time}_s{ssize}_r{rec}_{rep}.trees'),
		demes_file = 'results/simulations/scenario_{sc}.json',
	output:
		trees_file = 'results/simulations/sim_slim_sel_rep/sim_slim_sel_inter_{sc}_{type}_t{time}_s{ssize}_r{rec}_{rep}.trees',
	params:
		neutral_mut_rate = 1e-08,
	resources:
		mem_mb = 3_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/sim_slim_postprocessing.py'