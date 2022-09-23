import numpy as np
import tskit
import matplotlib.pyplot as plt
import seaborn as sns
import dabest
import pandas as pd


def ts_freq(ts, sample_sets=None):
	if sample_sets is None:
		sample_sets = [ts.samples()]
	n = np.array([len(x) for x in sample_sets])

	def f(x):
		return x / n

	return ts.sample_count_stat(
		sample_sets,
		f,
		len(sample_sets),
		span_normalise=False,
		windows="sites",
		polarised=True,
		mode="site",
		strict=False,
	).T


def ts_delete_recurrent_mutations(ts):
	if ts.num_mutations == 0:
		raise ValueError("treesequence does not have mutations, nothing to do")
	tables = ts.dump_tables()
	tables.mutations.clear()
	for site in ts.sites():
		tables.mutations.append(site.mutations[0])
		# because they are ordered by time, keep oldest
	return tables.tree_sequence()


# def variant_pseudohaploid_freq(variant, base_indices, N, rng):
# 	return np.mean(variant.genotypes[base_indices + rng.integers(0, 2, N)] > 0)

# def ts_freq_pseudohaploid_draw(ts, inds, rng):
# 	N = len(inds)
# 	base_indices = np.array(range(0, N*2, 2))
# 	base_nodes = np.array([ts.individual(i).nodes for i in inds]).flatten()
# 	freqs = [
# 		variant_pseudohaploid_freq(var, base_indices, N, rng) for var in ts.variants(samples=base_nodes)
# 	]
# 	return np.array(freqs)


def freq_pseudohaploid_draw(gen_mat, rng):
	# assumes gen_mat is ordered such that each node of same ind
	# are adjacent on axis 1
	N = gen_mat.shape[1]
	base_indices = np.array(range(0, N, 2))
	indices = rng.integers(0, 2, size=(gen_mat.shape[0], int(N / 2))) + base_indices
	gen_subset = np.take_along_axis(gen_mat, indices, axis=1)
	return np.mean(gen_subset > 0, axis=1)


def geno_pseudohaploid_draw(gen_mat, rng):
	# assumes gen_mat is ordered such that each node of same ind
	# are adjacent on axis 1
	N = gen_mat.shape[1]
	base_indices = np.array(range(0, N, 2))
	indices = rng.integers(0, 2, size=(gen_mat.shape[0], int(N / 2))) + base_indices
	gen_subset = np.take_along_axis(gen_mat, indices, axis=1)
	return gen_subset


def ts_freq_wrap_sets(ts, ind_sets: list, rng):
	af = np.zeros((len(ind_sets), ts.num_sites))
	gen_mat = ts.genotype_matrix()
	for i, x in enumerate(ind_sets):
		set_nodes = np.array([ts.individual(i).nodes for i in x]).flatten()
		af[i] = freq_pseudohaploid_draw(gen_mat[:, set_nodes], rng)
	return af


def ts_freq_wrap_samples_sets(ts, samples_sets: list, rng):
	af = np.zeros((len(samples_sets), ts.num_sites))
	gen_mat = ts.genotype_matrix()
	for i, x in enumerate(samples_sets):
		af[i] = freq_pseudohaploid_draw(gen_mat[:, x], rng)
	return af


# def ts_freq_wrap_sets_thinning(ts, ind_sets, to_drop, rng):
#     # to_drop = np.where(rng.uniform(size=ts.num_sites) > prop)[0]
#     thinned_ts = ts.delete_sites(to_drop)
#     gen_mat = thinned_ts.genotype_matrix()
#     af = np.zeros((len(ind_sets), thinned_ts.num_sites))
#     for i, x in enumerate(ind_sets):
#         set_nodes = np.array([ts.individual(i).nodes for i in x]).flatten()
#         af[i] = ts_freq_pseudohaploid_draw(gen_mat[:, set_nodes], rng)
#     return af


def ts_freq_wrap_sets_2(ts, times, n_sample, rng):
	af = np.zeros((len(times), ts.num_sites))
	gen_mat = ts.genotype_matrix()
	nodetable = ts.tables.nodes
	for i, t in enumerate(times):
		set_nodes = rng.choice(
			np.where((nodetable.time == t) & (nodetable.population == 3))[0],
			n_sample,
			replace=False,
		)
		af[i] = np.mean(gen_mat[:, set_nodes] > 0, axis=1)
	return af


def ts_individuals_at(ts, time, population):
	has_indiv = ts.tables.nodes.individual >= 0
	which_indiv = ts.tables.nodes.individual[has_indiv]
	individual_times = np.zeros(ts.num_individuals)
	individual_times[which_indiv] = ts.tables.nodes.time[has_indiv]
	individual_populations = np.repeat(np.int32(-1), ts.num_individuals)
	individual_populations[which_indiv] = ts.tables.nodes.population[has_indiv]
	res_inds = np.where(
		(individual_times == time) & (individual_populations == population)
	)[0]
	return res_inds


def ts_get_times(ts):
	has_indiv = ts.tables.nodes.individual >= 0
	which_indiv = ts.tables.nodes.individual[has_indiv]
	individual_times = np.zeros(ts.num_individuals)
	individual_times[which_indiv] = ts.tables.nodes.time[has_indiv]
	return np.unique(individual_times)


def ts_draw_ind_sets(ts, times, rng, pop, n_sample):
	ind_sets = [
		rng.choice(
			ts_individuals_at(ts, t, pop),
			size=n_sample,
			replace=False,
		)
		for t in times
	]
	return ind_sets


def ts_draw_samples_sets(ts, times: list, rng, pop, n_sample):
	nodes = ts.tables.nodes
	all_samples = [
		(nodes.population == pop) & (nodes.time == t) & (nodes.flags == 1)
		for t in times
	]
	ind_id = [np.unique(nodes[s].individual) for s in all_samples]
	# take only n_sample inds
	# take both nodes of each ind
	samples_sets = [
		np.where(
			s & np.isin(nodes.individual, rng.choice(id, size=n_sample, replace=False))
		)[0]
		for s, id in zip(all_samples, ind_id)
	]
	return samples_sets


def ts_draw_ind_sets_ref(ts, time, rng, pops, n_sample):
	ind_sets = [
		rng.choice(
			ts_individuals_at(ts, time, pop),
			size=n_sample,
			replace=False,
		)
		for pop in pops
	]
	return ind_sets


def get_list_freqs(list_ts, list_ind_sets, rng):
	af = []
	for ts, x in zip(list_ts, list_ind_sets):
		af.append(ts_freq_wrap_sets(ts, x, rng))
	return af


def compute_deltas(af):
	# assume af is in increasing time order
	# flipping ensures we do x_{t+1} - x_{t}
	deltas = np.flip(
		np.diff(np.flip(af, axis=0), n=1, axis=0),
		axis=0,
	)
	return deltas


def compute_cov(deltas):
	cov = np.cov(deltas, bias=True)
	return cov


def kth_diag_indices(a, k=0):
	rows, cols = np.diag_indices_from(a)
	if k < 0:
		return rows[-k:], cols[:k]
	elif k > 0:
		return rows[:-k], cols[k:]
	else:
		return rows, cols


def compute_bias_vector(af, n_sample):
	hh = af * (1 - af)
	bias = np.mean(1 / (n_sample - 1) * hh, axis=1)
	return bias


def create_bias_correction_matrix(b):
	size = len(b) - 1
	c = np.zeros((size, size))
	c[kth_diag_indices(c, k=0)] += b[:-1] + b[1:]
	c[kth_diag_indices(c, k=1)] -= b[1:-1]
	c[kth_diag_indices(c, k=-1)] -= b[1:-1]
	return c


def compute_bias_covmat(af, n_sample):
	b = compute_bias_vector(af, n_sample)
	c = create_bias_correction_matrix(b)
	return c


def apply_sampling_correction(a, B):
	c = a.copy()
	c[kth_diag_indices(c, k=0)] -= B[:-1] + B[1:]
	c[kth_diag_indices(c, k=1)] += B[1:-1]
	c[kth_diag_indices(c, k=-1)] += B[1:-1]
	return c


# =================================================
# Admixtures

# def ts_get_admixture_proportions_nodes(ts, admix_nodes, popA_nodes, popB_nodes, popC_nodes):
# 	admix_proportions = np.zeros((len(admix_nodes), 3))
# 	edges = ts.tables.link_ancestors(admix_nodes, np.concatenate((popA_nodes, popB_nodes, popC_nodes)))
# 	for ix, ind in enumerate(admix_nodes):
# 		A_edges = edges[np.isin(edges.child, admix_nodes) & np.isin(edges.parent, popA_nodes)]
# 		B_edges = edges[np.isin(edges.child, admix_nodes) & np.isin(edges.parent, popB_nodes)]
# 		C_edges = edges[np.isin(edges.child, admix_nodes) & np.isin(edges.parent, popC_nodes)]
# 		span_A = np.sum(A_edges.right - A_edges.left)
# 		span_B = np.sum(B_edges.right - B_edges.left)
# 		span_C = np.sum(C_edges.right - C_edges.left)
# 		admix_proportions[ix,0] = span_A/(span_A + span_B + span_C)
# 		admix_proportions[ix,1] = span_B/(span_A + span_B + span_C)
# 		admix_proportions[ix,2] = span_C/(span_A + span_B + span_C)
# 	return admix_proportions


def ts_get_admixture_proportions_inds(
	ts, admix_inds, ancestral_pops_nodes: list[np.array]
):
	N_anc_pop = len(ancestral_pops_nodes)
	admix_proportions = np.zeros((len(admix_inds), N_anc_pop))
	admix_nodes = np.array([ts.individual(ind).nodes for ind in admix_inds]).flatten()
	edges = ts.tables.link_ancestors(admix_nodes, np.concatenate(ancestral_pops_nodes))
	for ix, ind in enumerate(admix_inds):
		child_nodes = ts.individual(ind).nodes
		anc_edges = [
			edges[np.isin(edges.child, child_nodes) & np.isin(edges.parent, pop_nodes)]
			for pop_nodes in ancestral_pops_nodes
		]
		spans = [np.sum(pop_edges.right - pop_edges.left) for pop_edges in anc_edges]
		for j in range(N_anc_pop):
			admix_proportions[ix, j] = spans[j] / np.sum(spans)
	return admix_proportions


def ts_get_admixtures(ts, sample_sets, census_time, ancestral_pops):
	anc_pop_nodes = [
		np.where(
			(ts.tables.nodes.population == p) & (ts.tables.nodes.time == census_time)
		)[0]
		for p in ancestral_pops
	]
	admixtures = []
	for g in sample_sets:
		admixtures.append(ts_get_admixture_proportions_inds(ts, g, anc_pop_nodes))
	return np.stack(admixtures)


def compute_exp_admix_freq(admixtures, ref_af):
	admix_af = np.dot(np.mean(admixtures, axis=1), ref_af)
	return admix_af


# ==================================================
# Overall function


def tree_analysis(
	ts: tskit.TreeSequence,
	times: list,
	census_time,
	n_sample=20,
	samples_sets: list[np.array] = None,
	admix_pop: int = 3,
	ref_pops: list[int] = [1, 2],
	rng=None,
):
	if rng is None:
		rng = np.random.default_rng()

	if samples_sets is None:
		samples_sets = ts_draw_samples_sets(ts, times, rng, admix_pop, n_sample)

	af = ts_freq_wrap_samples_sets(ts, samples_sets, rng)
	deltas = compute_deltas(af)
	cov = compute_cov(deltas)
	bias = compute_bias_covmat(af, n_sample)

	ind_sets = [np.unique([ts.node(u).individual for u in s]) for s in samples_sets]
	admix = ts_get_admixtures(ts, ind_sets, census_time, ref_pops)

	ref_nodes = [
		np.where(
			(ts.tables.nodes.population == p)
			& (ts.tables.nodes.time == census_time)
			& (ts.tables.nodes.flags == 1)
		)[
			0
		]  # taking all available samples
		for p in ref_pops
	]

	ref_af = ts_freq(ts, ref_nodes)

	admix_af = compute_exp_admix_freq(admix, ref_af)
	admix_deltas = compute_deltas(admix_af)
	admix_cov = compute_cov(admix_deltas)

	return {
		"samples_sets": samples_sets,
		"af": af,
		"cov": cov,
		"bias": bias,
		"admix": admix,
		"admix_cov": admix_cov,
		"admix_af": admix_af,
	}


# =================================================
# Plotting


def dabest_bootstrap_test(control: np.array, test: np.array):
	# returns True if 0 is not in 95% CI, False otherwise
	dim = control.shape[-1]
	res = np.zeros((dim, dim), dtype=bool)
	for i in range(dim):
		for j in range(i + 1):
			dab = dabest.load(
				data=pd.DataFrame({"control": control[:, i, j], "test": test[:, i, j]}),
				idx=("control", "test"),
				resamples=5000,
			)
			t = (
				dab.mean_diff.statistical_tests.bca_low[0]
				* dab.mean_diff.statistical_tests.bca_high[0]
				>= 0
			)
			if t:
				res[i, j] = True
	res = res | res.T  # cp lower triangle to upper
	return res


def plot_covmats(
	list_covmat,
	list_pval=None,
	list_titles=None,
	main_title=None,
	delta_labels=None,
	scales=None,
):
	N = len(list_covmat)
	N_delta = list_covmat[0].shape[0]
	mask_covmat = np.arange(N_delta)[:, None] > np.arange(N_delta)

	if delta_labels is None:
		delta_labels = [f"$\Delta_{{{x}}}$" for x in range(N_delta)]

	fig, axs = plt.subplots(1, N, figsize=(5 * N + 5, 5))
	for i in range(N):
		tmp_mat = list_covmat[i].copy()
		tmp_mat[mask_covmat] = np.nan
		if scales is None:
			scale_max = np.max(np.abs([np.nanmin(tmp_mat), np.nanmax(tmp_mat)]))
		else:
			scale_max = scales[i]
		sns.heatmap(
			tmp_mat.T,
			cmap="vlag",
			vmin=-scale_max,
			vmax=scale_max,
			xticklabels=delta_labels,
			yticklabels=delta_labels,
			linewidths=0.5,
			ax=axs[i],
		)
		if list_titles is not None:
			axs[i].set_title(list_titles[i])

		if list_pval is not None:
			if list_pval[i] is not None:
				sig = list_pval[i]
				if sig.dtype == float:
					Bonf_alpha = 0.05 / (N_delta * (N_delta - 1) / 2 + N_delta)
					sig = list_pval[i]
					for j in range(sig.shape[0]):
						for z in range(j, sig.shape[0]):
							if sig[j, z] < Bonf_alpha:
								_ = axs[i].text(
									j + 0.5, z + 0.5, "*", ha="center", va="center"
								)
				if sig.dtype == bool:
					for j in range(sig.shape[0]):
						for z in range(j, sig.shape[0]):
							if sig[j, z]:
								_ = axs[i].text(
									j + 0.5, z + 0.5, "*", ha="center", va="center"
								)

	fig.subplots_adjust(bottom=0.2)
	if main_title is not None:
		fig.suptitle(main_title)
	return fig


# =======================================================
# Computing G


def compute_G(af, cov, bias, admix_cov, n_sample):
	af1 = af[-1]
	af0 = af[0]
	b = np.mean(1/(n_sample - 1)*af1*(1-af1)) + np.mean(1/(n_sample - 1)*af0*(1-af0))
	tot_var = np.var(af1 - af0) - b
	C = cov - bias - admix_cov
	num = np.sum(C - np.diag(np.diag(C)))  # remove the variances
	G = num / tot_var
	return G


def compute_num_G(C):
	num = np.sum(C - np.diag(np.diag(C)))
	return num


def compute_tot_var(af, n_sample):
	af1 = af[-1]
	af0 = af[0]
	b = np.mean(1/(n_sample - 1)*af1*(1-af1)) + np.mean(1/(n_sample - 1)*af0*(1-af0))
	tot_var = np.var(af1 - af0) - b
	return tot_var


def run_bootstrap(tiled_stat, N_bootstrap, weights):
	straps = []
	for b in np.arange(N_bootstrap):
		bidx = np.random.randint(0, len(tiled_stat), size=len(tiled_stat))
		straps.append(
			np.average(np.stack(tiled_stat)[bidx], axis=0, weights=weights[bidx])
		)
	return np.stack(straps)


def run_bootstrap_ratio(tiled_stat_num, tiled_stat_denom, N_bootstrap, weights):
	straps_num = []
	straps_denom = []
	for b in np.arange(N_bootstrap):
		bidx = np.random.randint(0, len(tiled_stat_num), size=len(tiled_stat_num))
		straps_num.append(
			np.average(np.stack(tiled_stat_num)[bidx], axis=0, weights=weights[bidx])
		)
		straps_denom.append(
			np.average(np.stack(tiled_stat_denom)[bidx], axis=0, weights=weights[bidx])
		)
	return np.stack(straps_num)/np.stack(straps_denom)


def bootstrap(
	ts,
	af,
	admix_af,
	n_sample,
	N_bootstrap=5e3,
	tile_size: int = int(1e6),
):

	tiles = [(i, i + tile_size) for i in range(0, int(ts.sequence_length), tile_size)]
	sites = ts.tables.sites
	tile_masks = [
		(start <= sites.position) & (sites.position < stop) for start, stop in tiles
	]
	n_loci = np.array([np.sum(tile) for tile in tile_masks])

	deltas = compute_deltas(af)
	admix_deltas = compute_deltas(admix_af)

	tiled_af = [af[:, mask] for mask in tile_masks]

	tiled_cov = [compute_cov(deltas[:, mask]) for mask in tile_masks]

	tiled_bias = [compute_bias_covmat(a, n_sample) for a in tiled_af]

	tiled_admix_cov = [compute_cov(admix_deltas[:, mask]) for mask in tile_masks]

	tiled_corr_cov = [
		c - b - a
		for c, b, a in zip(tiled_cov, tiled_bias, tiled_admix_cov)
	]

	tiled_num_G = [
		compute_num_G(c - b - ac)
		for c, b, ac in zip(tiled_cov, tiled_bias, tiled_admix_cov)
	]
	tiled_tot_var = [compute_tot_var(a, n_sample) for a, b in zip(tiled_af, tiled_bias)]

	tiled_num_Ap = [np.sum(ac) for ac in tiled_admix_cov]

	weights = n_loci / np.sum(n_loci)
	straps_cov = run_bootstrap(tiled_cov, N_bootstrap, weights)
	straps_bias = run_bootstrap(
		[c - b for c, b in zip(tiled_cov, tiled_bias)], 
		N_bootstrap, weights
	)
	straps_corr_cov = run_bootstrap(tiled_corr_cov, N_bootstrap, weights)
	straps_G = run_bootstrap_ratio(tiled_num_G, tiled_tot_var, N_bootstrap, weights)
	straps_Ap = run_bootstrap_ratio(tiled_num_Ap, tiled_tot_var, N_bootstrap, weights)

	return (straps_cov, straps_bias, straps_corr_cov, straps_G, straps_Ap)
