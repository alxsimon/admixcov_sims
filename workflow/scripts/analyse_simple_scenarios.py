#%%
import admixcov as ac
import tskit
import demes
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

#%%
trees_file = snakemake.input['trees_file']
demes_file = snakemake.input['demes_file']
unit_n_sample = snakemake.params['n_sample']
unit_ref_n_sample = snakemake.params['ref_n_sample']
census_time = snakemake.params['census_time']
info = snakemake.output['info']

drop_times = 2 if 'slim' in trees_file else 1

ts = tskit.load(trees_file)
times = np.flip(ac.ts.get_times(ts))[drop_times:]
n_samples = [unit_n_sample]*len(times)
graph = demes.load(demes_file)
admix_pop = ts.num_populations - 1
refs = [
    {'pop': i, 'time': census_time, 'n': unit_ref_n_sample}
    for i in range(admix_pop)
]
rng = np.random.default_rng()

with open(info, 'w') as f:
    print(ts, file=f)
    print('\n', file=f)
    print('Analysed times', file=f)
    print(times, file=f)
    print('\n', file=f)


samples_nodes = ac.ts.draw_sample_sets(ts, times, rng, admix_pop, n_samples)
ref_nodes = [
    ac.ts.draw_sample_set_at(ts, r['time'], rng, r['pop'], r['n'])
    for r in refs
]

ref_geno_sets, geno_sets = ac.ts.get_pseudohap_genotypes(
    ts, rng,
    samples=(ref_nodes, samples_nodes),
)

af = np.stack([
    np.mean(G, axis=1)
    for G in geno_sets
], axis=0)
ref_af = np.stack([
    np.mean(G, axis=1)
    for G in ref_geno_sets
], axis=0)

# admixture computations
samples_inds = [np.unique([ts.node(u).individual for u in s]) for s in samples_nodes]
ancestral_census_nodes = [
    np.where(
        (ts.tables.nodes.population == r['pop'])
        & (ts.tables.nodes.time == r['time'])
    )[0]
    for r in refs
]
Q = np.stack([
    np.mean(
        ac.ts.get_admixture_proportions(ts, g, ancestral_census_nodes),
        axis=0,
    )
    for g in samples_inds
])

alpha_mask = np.array(
    [p.proportions for p in graph.pulses]
) > 0 # create alphas from graph
alphas = ac.q2a_simple(Q, alpha_mask)
alphas[alphas < 0] = 0 # issues due to drift
covmat = ac.get_covariance_matrix(
    af, bias=True,
    sample_size=np.array(n_samples),
)
admix_cov = ac.get_admix_covariance_matrix(
    Q, ref_af, bias=True,
    ref_sample_size=np.array([d['n'] for d in refs]),
)
var_drift = ac.solve_for_variances(
    np.diag(covmat - admix_cov), alphas,
)
drift_err = ac.get_drift_err_matrix(
    var_drift, alphas,
)


#%%
k = covmat.shape[0]
G = []
G_nc = []
Ap = []
totvar = []
totvar_adj = []
for i in range(1, k + 1):
    totvar.append(np.sum(covmat[:i, :i]))
    totvar_adj.append(np.sum((covmat - admix_cov - drift_err)[:i, :i]))
    G.append(
        ac.get_matrix_sum(
            (covmat - admix_cov - drift_err)[:i, :i],
            include_diag=False, abs=False
        ) / totvar[-1]
    )
    G_nc.append(
        ac.get_matrix_sum(
            covmat[:i, :i],
            include_diag=False, abs=False
        ) / totvar[-1]
    )
    Ap.append(
        ac.get_matrix_sum(
            admix_cov[:i, :i],
            include_diag=True, abs=False
        ) / totvar[-1]
    )

with open(info, 'a') as f:
    print('Total variance:', file=f)
    print(totvar, file=f)
    print('G:', file=f)
    print(G, file=f)
    print('G non-corrected:', file=f)
    print(G_nc, file=f)
    print('A\':', file=f)
    print(Ap, file=f)
    print('Total variance adjusted:', file=f)
    print(totvar_adj, file=f)
    print('\n', file=f)


#%%
tile_idxs = ac.ts.create_tile_idxs(ts, tile_size=int(5e5))

n_sample = np.array([n_samples]*af.shape[1]).T
n_sample_ref = np.array([[d['n'] for d in refs]]*af.shape[1]).T
tiled_af = [af[:, idx] for idx in tile_idxs]
tiled_sample_size = [n_sample[:, idx] for idx in tile_idxs]

assert af.shape == n_sample.shape
tiled_cov = np.stack([
    ac.get_covariance_matrix(a, bias=True, sample_size=n)
    for a, n in zip(tiled_af, tiled_sample_size)
])

assert ref_af.shape == n_sample_ref.shape
tiled_admix_cov = np.stack([
    ac.get_admix_covariance_matrix(
        Q, ref_af[:, idx], bias=True,
        ref_sample_size=n_sample_ref[:, idx],
    )
    for idx in tile_idxs
])

tiled_drift_err = np.stack([
    ac.get_drift_err_matrix(
        ac.solve_for_variances(np.diag(c - a), alphas),
        alphas,
    )
    for c, a in zip(tiled_cov, tiled_admix_cov)
])

tiled_corr_cov = np.stack([
    c - a - d for c, a, d in zip(tiled_cov, tiled_admix_cov, tiled_drift_err)
])

n_loci = np.array([tile.size for tile in tile_idxs])
weights = n_loci / np.sum(n_loci)

# do the bootstraps
straps_cov_nc = ac.bootstrap_stat(tiled_cov, weights, 1e4)
straps_cov = ac.bootstrap_stat(tiled_corr_cov, weights, 1e4)

straps_G = []
straps_G_nc = []
straps_Ap = []
k = tiled_cov.shape[1]
for i in range(1, k + 1):
    tmp_totvar = np.sum(tiled_cov[:, :i, :i], axis=(1, 2))
    straps_G.append(
        ac.bootstrap_ratio(
            np.stack([ac.get_matrix_sum(c) for c in tiled_corr_cov[:, :i, :i]]),
            tmp_totvar,
            weights,
            1e4,
            statistic=G[i - 1],
        )
    )
    straps_G_nc.append(
        ac.bootstrap_ratio(
            np.stack([ac.get_matrix_sum(c) for c in tiled_cov[:, :i, :i]]),
            tmp_totvar,
            weights,
            1e4,
            statistic=G_nc[i - 1],
        )
    )
    straps_Ap.append(
        ac.bootstrap_ratio(
            np.stack([ac.get_matrix_sum(c, include_diag=True) for c in tiled_admix_cov[:, :i, :i]]),
            tmp_totvar,
            weights,
            1e4,
            statistic=Ap[i - 1],
        )
    )


#%%
d = {
    'times': times,
    'Q': Q,
    'covmat': covmat,
    'admix_cov': admix_cov,
    'drift_err': drift_err,
    'G': G,
    'Ap': Ap,
}


# %%
time_padding = 10

def plot_ci_line(x, CI, ax, color='black', linestyle='-', marker='o'):
    lower = CI[0]
    m =CI[1]
    upper = CI[2]
    yerr = np.array([m - lower, upper - m])
    ax.errorbar(x, m, yerr, color=color, ecolor=color, marker=marker, linestyle=linestyle)

def cov_lineplot(times, CIs, ax, colors, linestyle='-', marker='o'):
    k = len(times) - 1
    for i in range(k-1):
        plot_ci_line(np.array(times[i+1:-1]) + 1, np.stack(CIs)[:, i, i+1:], ax, color=colors[i], linestyle=linestyle, marker=marker)
    ax.hlines(y=0, xmin=0, xmax=times[1] + time_padding, linestyles='dotted', colors='black')
    ax.set_xlim(times[1] + time_padding, -time_padding)
    ax.set_xlabel('time')
    ax.set_ylabel('covariance')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def combine_covmat_CIs(ci_l, ci_u):
	N = ci_l[1].shape[0]
	res = tuple([x.copy() for x in ci_l])
	tri_up = np.triu_indices(N, k=1)
	for k in range(3):
		res[k][np.arange(N)[:, None] == np.arange(N)] = np.nan
		res[k][tri_up] = ci_u[k][tri_up]
	return res

def plot_covmat(CI, ax, scale_max=None):
    N_delta = CI[1].shape[0]
    delta_labels = [f"$\\Delta_{{{x}}}$" for x in range(N_delta)]
    tmp_mat = CI[1].copy()
    if scale_max is None:
        scale_max = np.max(np.abs([np.nanmin(tmp_mat), np.nanmax(tmp_mat)]))
    sns.heatmap(
        tmp_mat.T,
        cmap="vlag",
        vmin=-scale_max,
        vmax=scale_max,
        xticklabels=delta_labels,  # type: ignore
        yticklabels=delta_labels,  # type: ignore
        linewidths=0.5,  # type: ignore
        ax=ax,
    )
    sig = (CI[0] * CI[2]) > 0
    for z in range(sig.shape[0]):
        for j in range(sig.shape[0]):
            if (sig[j, z]) & (z != j):
                _ = ax.text(
                    j + 0.5, z + 0.5, "*", ha="center", va="center"
                )
    ax.axline(xy1=(N_delta, N_delta), slope=1, color='black')


colors_oi = [
    '#000000', # black
    '#D55E00', # vermillion
    '#0072B2', # blue
    '#009E73', # green
    '#E69F00', # orange
    '#56B4E9', # sky blue
    '#CC79A7', # pink
    '#F0E442', # yellow
]

fig, axs = plt.subplots(3, 2, figsize=(10, 8))

for i in range(d['Q'].shape[1]):
	axs[0,0].plot(d['times'], d['Q'][:, i], '-o', label=f"Pop{i}")
axs[0,0].set_xlim(d['times'][0] + time_padding, -time_padding)
axs[0,0].set_ylim((0,1))
axs[0,0].set_ylabel("Mean ancestry")
axs[0,0].set_xlabel("Time point")
axs[0,0].legend(loc="upper right")

combined_ci = combine_covmat_CIs(straps_cov, straps_cov_nc)
scale_max = np.max(np.abs([np.nanmin(combined_ci[1] - np.diag(np.diag(combined_ci[1]))), np.nanmax(combined_ci[1] - np.diag(np.diag(combined_ci[1])))]))
plot_covmat(combined_ci, axs[0, 1], scale_max)
axs[0,1].set_title('covariance matrix (raw lower, corrected upper)')

ymin, ymax = (
    np.min(np.concatenate([np.tril(straps_cov_nc[0]), np.tril(straps_cov[0])])),
    np.max(np.concatenate([np.tril(straps_cov_nc[2]), np.tril(straps_cov[2])])),
)
cov_lineplot(d['times'], straps_cov_nc, axs[1, 0], colors=colors_oi)
axs[1, 0].set_ylim(ymin, ymax)
axs[1, 0].set_ylabel('raw covariance (without bias)')
cov_lineplot(d['times'], straps_cov, axs[1, 1], colors=colors_oi)
axs[1, 1].set_ylim(ymin, ymax)
axs[1, 1].set_ylabel('admixture corrected covariance')

sns.lineplot(x=d['times'][1:], y=totvar, ax=axs[2, 0], marker='o')
axs[2, 0].set_xlim(d['times'][1] + time_padding, -time_padding)
axs[2, 0].set_ylim(0)
axs[2, 0].set_ylabel('Total variance (t)')

plot_ci_line(d['times'][1:], np.stack(straps_G_nc).T, ax=axs[2, 1], linestyle='dashed')
plot_ci_line(d['times'][1:], np.stack(straps_G).T, ax=axs[2, 1])
plot_ci_line(d['times'][1:], np.stack(straps_Ap).T, ax=axs[2, 1], marker='s', color='blue')
ymin, ymax = (
    np.min(np.concatenate([np.stack(straps_G_nc), np.stack(straps_G), np.stack(straps_Ap)])),
    np.max(np.concatenate([np.stack(straps_G_nc), np.stack(straps_G), np.stack(straps_Ap)])),
)
axs[2, 1].set_xlim(d['times'][1] + time_padding, -time_padding)
axs[2, 1].set_ylim(ymin, ymax)
axs[2, 1].hlines(y=0, xmin=0, xmax=d['times'][1] + time_padding, linestyles='dotted', colors='black')
axs[2, 1].set_xlabel('time')
axs[2, 1].set_ylabel("G(t) or A'(t)")

fig.tight_layout()

fig.savefig(
    snakemake.output['main_fig'],
)

#%%
if 'slim' in trees_file:
    fig, ax = plt.subplots(figsize=(6, 4))
    ztb = pd.read_csv(trees_file.replace('.trees', '_pheno.tsv'), sep='\t')
    ztb['bgen'] = ztb.gen.max() - ztb.gen
    sns.lineplot(
        ztb[ztb.bgen < times[0] + 20], 
        x='bgen', y='mean_z', style='pop', hue='pop',
        ax=ax,
    )
    ax.set_xlim(xmin=times[0] + 20, xmax=-5)
    ax.set_xlabel('generations')
    ax.set_ylabel('mean phenotype')
    fig.savefig(snakemake.output['main_fig'].replace('.pdf', '_pheno.pdf'))