
rule analyse_msprime_simple_scenarios:
	input:
		files = expand(
			'results/simulations/sim_msprime_rep/sim_msprime_scenario_{{sc}}_{rep}.trees',
			rep=range(config['N_rep']),
		),
		demes_file = 'results/simulations/scenario_{sc}.yaml',
	output:
		pickle = 'results/simulations/sim_msprime_scenario_{sc}.pickle',
	params:
		census_time = 200,
		n_sample = 30,
		ref_n_sample = 30,
	resources:
		mem_mb = 10_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/analyse_simple_scenarios.py'


rule analyse_slim_sel_simple_scenarios:
	input:
		files = expand(
			'results/simulations/sim_slim_sel_rep/sim_slim_sel_scenario_{{sc}}_{{type}}_t{{time}}_s{{ssize}}_{rep}.trees',
			rep=range(config['N_rep']),
		),
		demes_file = 'results/simulations/scenario_{sc}.yaml',
	output:
		pickle = 'results/simulations/sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}.pickle',
	params:
		census_time = 201,
		n_sample = 30,
		ref_n_sample = 30,
	resources:
		mem_mb = 10_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/analyse_simple_scenarios.py'


rule plot_msprime_simple_scenarios:
	input:
		pickle = 'results/simulations/sim_msprime_scenario_{sc}.pickle',
	output:
		main_fig = 'results/figures/sim_msprime_scenario_{sc}.pdf',
	resources:
		mem_mb = 4000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/plot_simple_scenarios.py'


rule plot_slim_sel_simple_scenarios:
	input:
		pickle = 'results/simulations/sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}.pickle',
	output:
		main_fig = 'results/figures/sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}.pdf',
		pheno_fig = 'results/figures/sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}_pheno.pdf'
	resources:
		mem_mb = 4000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/plot_simple_scenarios.py'


# for Patterson
# [37  69  26  23 273  38  62]
# [18, 21, 18]

rule analyse_msprime_europe_Patterson2022:
	input:
		files = expand(
			'results/simulations/sim_msprime_europe_Patterson2022/sim_msprime_europe_Patterson2022_{rep}.trees',
			rep=range(config['N_rep']),
		),
		demes_file = 'resources/model_europe_Patterson2022.yaml',
	output:
		pickle = 'results/simulations/sim_msprime_europe_Patterson2022.pickle',
	params:
		census_time = 160,
		n_samples = [37, 69, 26, 23, 273, 38, 62],
		ref_n_samples = [18, 21, 18],
	resources:
		mem_mb = 10_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/analyse_europe_Patterson2022.py'


rule plot_msprime_europe_Patterson2022:
	input:
		pickle = 'results/simulations/sim_msprime_europe_Patterson2022.pickle',
	output:
		main_fig = 'results/figures/sim_msprime_europe_Patterson2022.pdf',
	resources:
		mem_mb = 4000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/plot_europe_Patterson2022.py'