#!python
# import demes
# import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

try:
	snakemake
except:
	import workflow.scripts.funcs as fn
	print("Not running in Snakemake")
	print("Defining variables manually")
	# inputs
	trees_file = 'results/simulations/sim_eu_A.trees'
	rate_map_pickle = 'results/simulations/sim_eu_A_rate_map.pickle'
	# outputs
	fig_admix = 'results/tests/fig_admix.png'
	fig_covmat = 'results/tests/fig_covmat.png'
	fig_G_boot = 'results/tests/fig_G_boot.png'
	fig_A_boot = 'results/tests/fig_A_boot.png'
	fig_binned_rate = 'results/tests/fig_binned_rate.png'
	fig_PCA = 'results/tests/PCA.png'
	# params
	census_time = 210
	n_sample = 50
	bootstrap_tile_size = 5e5
else:
	import funcs as fn
	# inputs
	# demes_file = snakemake.input['demes_file']
	trees_file = snakemake.input['trees_file']
	rate_map_pickle = snakemake.input['rate_map_pickle']
	# outputs
	fig_admix = snakemake.output['fig_admix']
	fig_covmat = snakemake.output['fig_covmat']
	fig_G_boot = snakemake.output['fig_G_boot']
	fig_A_boot = snakemake.output['fig_A_boot']
	# fig_binned_rate = snakemake.output['fig_binned_rate']
	fig_PCA = snakemake.output['fig_PCA']
	# params
	census_time = snakemake.params['census_time']
	n_sample = snakemake.params['n_sample']
	bootstrap_tile_size = snakemake.params['bootstrap_tile_size']


ts = tskit.load(trees_file)

# Filter on MAF
inds_of_interest = np.concatenate([
	fn.ts_individuals_at(ts, 0, 0),
	# fn.ts_individuals_at(ts, 210, 0),
	# fn.ts_individuals_at(ts, 175, 1),
	# fn.ts_individuals_at(ts, 210, 2),
])
af_filter = fn.ts_freq(
	ts,
	[np.concatenate(
		[ts.individual(i).nodes for i in inds_of_interest]
	)]
).flatten()
dropped = np.where(af_filter < 0.05)[0]
ts = ts.delete_sites(dropped)

del af_filter, inds_of_interest, dropped
print(ts)

times = np.flip(fn.ts_get_times(ts))[1:]
# leave out first sample of ref
print(times)

# ==========================================
# Sampling nodes
rng = np.random.default_rng(66238645)
samples_nodes = fn.ts_draw_samples_sets(ts, times, rng, 0, n_sample)
refs = [
	{'pop': 0, 'time': 210},  # EEF
	{'pop': 1, 'time': 175},  # Yamnaya
	{'pop': 2, 'time': 210},  # WHG
]
ref_nodes = [
	np.where(
		(ts.tables.nodes.population == r['pop'])
		& (ts.tables.nodes.time == r['time'])
		& (ts.tables.nodes.flags == 1)
	)[
		0
	]  # taking all available samples
	for r in refs
]

# =========================================
# Getting pseudohaploid genotypes
concat_samples = np.concatenate(ref_nodes + samples_nodes)
geno = fn.ts_genotype_matrix_pseudohap(
	ts, rng, samples=concat_samples,
)

geno_sets = [
	geno[:, i:(i + n_sample)]
	for i in range(n_sample * 3, n_sample * 3 + n_sample * len(samples_nodes), n_sample)
]

ref_geno_sets = [
	geno[:, i:(i + n_sample)]
	for i in range(0, n_sample * 3, n_sample)
]


# =========================================
# admixture computations
samples_inds = [np.unique([ts.node(u).individual for u in s]) for s in samples_nodes]
ancestral_census_nodes = [
	np.where(
		(ts.tables.nodes.population == r['pop'])
		& (ts.tables.nodes.time == r['time'])
	)[0]
	for r in refs
]
admixtures = []
for g in samples_inds:
	admixtures.append(fn.ts_get_admixture_proportions_inds(ts, g, ancestral_census_nodes))
admixtures = np.stack(admixtures)

mean_admix = np.mean(admixtures, axis=1)
fig, ax = plt.subplots()
ax.plot(times, mean_admix[:, 0], "o-", label="Anatolian", color='orange')
ax.plot(times, mean_admix[:, 1], "^-", label="Yamnaya", color='green')
ax.plot(times, mean_admix[:, 2], "s-", label="WHG", color='red')
_ = ax.set_ylim((0, 1))
_ = ax.set_xlabel("Time ago (generations)")
_ = ax.set_ylabel("Mean ancestry")
ax.legend(loc='upper right')
plt.gca().invert_xaxis()
fig.savefig(fig_admix)

# =========================================
# Computing freqs and covariance
# af = fn.ts_freq_wrap_samples_sets(ts, samples_nodes, rng)
af = np.stack([
	np.mean(geno > 0, axis=1)
	for geno in geno_sets
], axis=0)
cov = fn.compute_cov(
	fn.compute_deltas(af)
)
bias = fn.compute_bias_covmat(af, n_sample)

ref_af = np.stack([
	np.mean(geno > 0, axis=1)
	for geno in ref_geno_sets
], axis=0)

admix_af = fn.compute_exp_admix_freq(admixtures, ref_af)
admix_cov = fn.compute_cov(
	fn.compute_deltas(admix_af)
)

G = fn.compute_G(af, cov, bias, admix_cov, n_sample)
print(G)
Ap = np.sum(admix_cov) / fn.compute_tot_var(af, n_sample)
print(Ap)

straps_cov, straps_bias, straps_corr_cov, straps_G, straps_Ap = fn.bootstrap(
	ts, af, admix_af, n_sample, N_bootstrap=1e4, tile_size=int(bootstrap_tile_size)
)

plt.figure()
_ = sns.histplot(straps_G)
plt.axvline(np.mean(straps_G), color='black')
plt.axvline(G, color='red')
plt.axvline(np.quantile(straps_G, q=0.025), linestyle='--')
plt.axvline(np.quantile(straps_G, q=0.975), linestyle='--')
plt.savefig(fig_G_boot)

plt.figure()
_ = sns.histplot(straps_Ap)
plt.axvline(np.mean(straps_Ap), color='black')
plt.axvline(Ap, color='red')
plt.axvline(np.quantile(straps_Ap, q=0.025), linestyle='--')
plt.axvline(np.quantile(straps_Ap, q=0.975), linestyle='--')
plt.savefig(fig_A_boot)


# ======================================
# PCA regression correction

N_loci, N_ind = geno.shape

p = ((1 + np.nansum(geno, axis=1)) / (2 + 2 * N_ind))
X = (geno.T - np.nanmean(geno, axis=1)) / ((p * (1 - p))**.5)
X[np.isnan(X)] = 0
del p

Psi = np.cov(X)
del X
e_values, e_vectors = np.linalg.eigh(Psi)
i = np.argsort(e_values)
i = i[::-1]
e_vectors = e_vectors[:, i]
e_values = e_values[i]

fig, ax = plt.subplots()
sns.scatterplot(
	x=e_vectors[(n_sample * 3):, 0],
	y=e_vectors[(n_sample * 3):, 1],
	ax=ax,
	label="EUR"
)
labels = ['EEF', 'Yam', 'WHG']
for i in range(3):
	sns.scatterplot(
		x=e_vectors[(n_sample * i):(n_sample * (i + 1)), 0],
		y=e_vectors[(n_sample * i):(n_sample * (i + 1)), 1],
		ax=ax,
		label=labels[i],
	)
fig.savefig(fig_PCA)


N_vec = 2
gamma = np.zeros((N_loci, N_vec))
for i in range(N_vec):
	gamma[:, i] = np.nansum(e_vectors[:, i] * geno, axis=1) / np.sum(e_vectors[:, i]**2)

geno_exp = np.dot(e_vectors[:, :N_vec], gamma.T)  # gamma_1*A1 + gamma_2*A2
geno_adj = geno.T - geno_exp
af_adj = np.stack([np.mean(geno_adj[i:(i + n_sample)], axis=0) for i in range(0, geno_adj.shape[0], n_sample)])
af_adj_ref = af_adj[:3]
af_adj = af_adj[3:]
cov_adj = fn.compute_cov(fn.compute_deltas(af_adj))

# ========================================
fig = fn.plot_covmats(
	[cov, cov - bias, cov - bias - admix_cov, cov_adj - bias, admix_cov],
	list_titles=[
		'raw covariances',
		'bias corrected',
		'bias and admixture corrected',
		'bias and admixture corrected (PCA)',
		'estimated admixture covariances',
	]
)
fig.savefig(fig_covmat)

# ==========================================
# Quantiles of recombination rates
# with open(rate_map_pickle, 'rb') as fr:
# 	rate_map = pickle.load(fr)

# rates = rate_map.rate[np.concatenate([
# 	np.where((rate_map.left <= s.position) & (rate_map.right > s.position))[0]
# 	for s in ts.sites()
# ])]

# q_r = np.quantile(rates, np.arange(0.1, 1, 0.1))
# rate_bins = np.digitize(rates, q_r)
# print(q_r)
# print(rate_bins.size)


# def compute_corrected_cov(af, admix_cov, n_sample):
# 	deltas = np.flip(
# 		np.diff(np.flip(af, axis=0), n=1, axis=0),
# 		axis=0,
# 	)
# 	covmat = np.cov(deltas)
# 	bias = fn.compute_bias_covmat(af, n_sample)
# 	return covmat - bias - admix_cov


# binned_cov_rate = []
# binned_G_rate = []
# binned_A_rate = []
# for i in range(9):
# 	bin_mask = (rate_bins == i)
# 	a_cov = np.cov(admix_deltas[:, bin_mask])
# 	binned_cov_rate.append(
# 		compute_corrected_cov(af[:, bin_mask], a_cov, n_sample)
# 	)

# 	tv = np.sum(binned_cov_rate[i])
# 	num = np.sum(binned_cov_rate[i] - np.diag(np.diag(binned_cov_rate[i])))
# 	binned_G_rate.append(num / tv)
# 	binned_A_rate.append(np.sum(a_cov) / tv)

# fig, axs = plt.subplots(2)
# axs[0].plot(range(9), binned_G_rate, 'o-')
# axs[0].set_title('G by rate bin')
# axs[1].plot(range(9), binned_A_rate, 'o-')
# axs[1].set_title("A' by rate bin")
# fig.tight_layout()
# fig.savefig(fig_binned_rate)

# ========================================
# Try estimating admixture estimation bias

# var_admix = np.var(admixtures, axis=1)
# mean_admix = np.mean(admixtures, axis=1)

# var_ref_noise = np.mean(1 / (n_sample - 1) * ref_af * (1 - ref_af), axis=1)

# admix_sampling_bias = np.sum((var_admix + mean_admix**2) * var_ref_noise, axis=1)
# admix_sampling_bias_covmat = fn.create_bias_correction_matrix(admix_sampling_bias)
