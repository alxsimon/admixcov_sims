
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
		mem_mb = 3_000,
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
		mem_mb = 3_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/analyse_simple_scenarios.py'


rule plot_msprime_simple_scenarios:
	input:
		pickle = 'results/simulations/sim_msprime_scenario_{sc}.pickle',
	output:
		main_fig = 'results/figures/fig_sim_msprime_scenario_{sc}.pdf',
	resources:
		mem_mb = 3_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/plot_simple_scenarios.py'


rule plot_slim_sel_simple_scenarios:
	input:
		pickle = 'results/simulations/sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}.pickle',
	output:
		main_fig = 'results/figures/fig_sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}.pdf',
		pheno_fig = 'results/figures/fig_sim_slim_sel_scenario_{sc}_{type}_t{time}_s{ssize}_pheno.pdf'
	resources:
		mem_mb = 3_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/plot_simple_scenarios.py'


# for Patterson
# [37  69  26  23 273  38  62]
# [18, 21, 18]

rule analyse_msprime_europe_uk:
	input:
		files = expand(
			'results/simulations/sim_msprime_europe_uk/sim_msprime_europe_uk_{rep}.trees',
			rep=range(config['N_rep']),
		),
		demes_file = 'resources/AncientEurope_4A21_mod.yaml',
	output:
		pickle = 'results/simulations/sim_msprime_europe_uk.pickle',
	params:
		census_time = 200,
		n_samples = [37, 69, 26, 23, 273, 38, 62],
		ref_n_samples = [18, 21, 18],
	resources:
		mem_mb = 9_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/analyse_europe_uk.py'


rule plot_msprime_europe_uk:
	input:
		pickle = 'results/simulations/sim_msprime_europe_uk.pickle',
		demes_file = 'resources/AncientEurope_4A21_mod.yaml',
	output:
		main_fig = 'results/figures/fig_sim_msprime_europe_uk.pdf',
		fig_demo = 'results/figures/fig_sim_msprime_europe_uk_demo.pdf',
	resources:
		mem_mb = 3_000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/plot_europe_uk.py'


def get_mem_mb_inter(wildcards):
	if wildcards.inter in ["20", "10", "5"]:
		mem = 6_000
	elif wildcards.inter in ["1"]:
		mem = 21_000
	else:
		mem = 3_000
	return mem

rule analyse_slim_sel_variable_interval:
	input:
		files = expand(
			'results/simulations/sim_slim_sel_rep/sim_slim_sel_inter_{{sc}}_{{type}}_t{{time}}_s{{ssize}}_i{{inter}}_{rep}.trees',
			rep=range(config['N_rep']),
		),
		demes_file = 'results/simulations/scenario_{sc}.yaml',
	output:
		pickle = 'results/simulations/sim_slim_sel_inter_{sc}_{type}_t{time}_s{ssize}_i{inter}.pickle',
	params:
		census_time = 201,
		n_sample = 30,
		ref_n_sample = 30,
	resources:
		mem_mb = get_mem_mb_inter,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/analyse_simple_scenarios.py'