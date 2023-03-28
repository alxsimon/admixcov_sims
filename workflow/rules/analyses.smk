
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
		mem_mb = 1000,
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
		mem_mb = 1000,
	conda:
		"../envs/popgensim.yaml"
	script:
		'../scripts/plot_simple_scenarios.py'


# for Patterson
# [37  69  26  23 273  38  62]
# [18, 21, 18]


# rule analyse_eu:
# 	input:
# 		trees_file = 'results/simulations/sim_eu_{variant}.trees',
# 		rate_map_pickle = 'results/simulations/sim_eu_{variant}_rate_map.pickle',
# 	output:
# 		fig_admix = 'results/figures/sim_eu_{variant}_admix.png',
# 		fig_covmat = 'results/figures/sim_eu_{variant}_covmat.png',
# 		fig_G_boot = 'results/figures/sim_eu_{variant}_G_boot.png',
# 		fig_A_boot = 'results/figures/sim_eu_{variant}_Ap_boot.png',
# 		# fig_binned_rate = 'results/figures/sim_eu_{variant}_binned_rate.png',
# 		fig_PCA = 'results/figures/sim_eu_{variant}_PCA.png',
# 	params:
# 		census_time = 210,
# 		n_sample = 50, # admix pop only
# 		bootstrap_tile_size = 5e5
# 	conda:
# 		"../envs/popgensim.yaml"
# 	script:
# 		'../scripts/sim_analysis_eu.py'


# rule analyse_simple_scenarios:
# 	input:
# 		trees_file = 'results/simulations/sim_scenario_{sc}.trees',
# 		demes_file = 'resources/scenario_{sc}.yaml',
# 	output:
# 		pickle = 'results/analyses/data_scenario_{sc}.pickle',
# 		text = 'results/analyses/stats_scenario_{sc}.txt',
# 	params:
# 		census_time = 200,
# 		n_sample = 20,
# 		ref_n_sample = 20,
# 	conda:
# 		"../envs/popgensim.yaml"
# 	script:
# 		'../scripts/analyse_simple_scenarios.py'

