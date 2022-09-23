#!python
# import demes
# import msprime
import tskit
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

try:
	snakemake
except:
	import workflow.scripts.funcs as fn
	print("Not running in Snakemake")
	print("Defining variables manually")
	# inputs
	# demes_file = '../../resources/model_europe.yaml'
	trees_file = 'results/simulations/europe_sim_noEHGpulse.trees'
	rate_map_pickle = 'results/simulations/europe_sim_rate_map_test.pickle'
	# outputs
	# params
	census_time = 210
	n_sample = 50
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
	fig_binned_rate = snakemake.output['fig_binned_rate']
	# params
	census_time = snakemake.params['census_time']
	n_sample = snakemake.params['n_sample']
	bootstrap_tile_size = snakemake.params['bootstrap_tile_size']


ts = tskit.load(trees_file)

times = np.flip(fn.ts_get_times(ts))[1:]
# leave out first sample of ref
print(times)


# for now analyze manually as I don't know how to get admixture with Yamnaya issue
rng = np.random.default_rng(276873265)
samples_sets = fn.ts_draw_samples_sets(ts, times, rng, 0, n_sample)
af = fn.ts_freq_wrap_samples_sets(ts, samples_sets, rng)
deltas = fn.compute_deltas(af)
cov = fn.compute_cov(deltas)
bias = fn.compute_bias_covmat(af, n_sample)

# ind_sets = [np.unique([ts.node(u).individual for u in s]) for s in samples_sets]
# admix = fn.ts_get_admixtures(ts, ind_sets, census_time, [1, 2, 3])

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
ref_af = fn.ts_freq(ts, ref_nodes)

# admixture
ind_sets = [np.unique([ts.node(u).individual for u in s]) for s in samples_sets]
anc_pop_nodes = [
	np.where(
		(ts.tables.nodes.population == r['pop'])
		& (ts.tables.nodes.time == r['time'])
	)[0]
	for r in refs
]
admixtures = []
for g in ind_sets:
	admixtures.append(fn.ts_get_admixture_proportions_inds(ts, g, anc_pop_nodes))
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

admix_af = fn.compute_exp_admix_freq(admixtures, ref_af)
admix_deltas = fn.compute_deltas(admix_af)
admix_cov = fn.compute_cov(admix_deltas)


fig = fn.plot_covmats(
	[cov, cov - bias, cov - bias - admix_cov, admix_cov],
	list_titles=[
		'raw covariances',
		'bias corrected',
		'bias and admixture corrected',
		'estimated admixture covariances',
	]
)
fig.savefig(fig_covmat)

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


# Quantiles of recombination rates
with open(rate_map_pickle, 'rb') as fr:
	rate_map = pickle.load(fr)

rates = rate_map.rate[np.concatenate([
	np.where((rate_map.left <= s.position) & (rate_map.right > s.position))[0]
	for s in ts.sites()
])]

q_r = np.quantile(rates, np.arange(0.1, 1, 0.1))
rate_bins = np.digitize(rates, q_r)
print(q_r)
print(rate_bins.size)


def compute_corrected_cov(af, admix_cov, n_sample):
	deltas = np.flip(
		np.diff(np.flip(af, axis=0), n=1, axis=0),
		axis=0,
	)
	covmat = np.cov(deltas)
	bias = fn.compute_bias_covmat(af, n_sample)
	return covmat - bias - admix_cov


binned_cov_rate = []
binned_G_rate = []
binned_A_rate = []
for i in range(9):
	bin_mask = (rate_bins == i)
	a_cov = np.cov(admix_deltas[:, bin_mask])
	binned_cov_rate.append(
		compute_corrected_cov(af[:, bin_mask], a_cov, n_sample)
	)

	tv = np.sum(binned_cov_rate[i])
	num = np.sum(binned_cov_rate[i] - np.diag(np.diag(binned_cov_rate[i])))
	binned_G_rate.append(num / tv)
	binned_A_rate.append(np.sum(a_cov) / tv)

fig, axs = plt.subplots(2)
axs[0].plot(range(9), binned_G_rate, 'o-')
axs[0].set_title('G by rate bin')
axs[1].plot(range(9), binned_A_rate, 'o-')
axs[1].set_title("A' by rate bin")
fig.tight_layout()
fig.savefig(fig_binned_rate)
