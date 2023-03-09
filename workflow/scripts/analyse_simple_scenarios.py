#%%
import admixcov as ac
import tskit
import demes
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

#%%
files = snakemake.input['files']
demes_file = snakemake.input['demes_file']
unit_n_sample = snakemake.params['n_sample']
unit_ref_n_sample = snakemake.params['ref_n_sample']
census_time = snakemake.params['census_time']
# info = snakemake.output['info']

drop_times = 2 if 'slim' in files[0] else 1

ts = tskit.load(files[0]) # extract info common to all trees
times = np.flip(ac.ts.get_times(ts))[drop_times:]
n_samples = [unit_n_sample]*len(times)
graph = demes.load(demes_file)
admix_pop = len(graph.demes) - 1
# not using ts.num_populations here as pyslim adds an additional one in ts
refs = [
    {'pop': i, 'time': census_time, 'n': unit_ref_n_sample}
    for i in range(admix_pop)
]
alpha_mask = np.array(
    [p.proportions for p in graph.pulses]
) > 0 # create alphas from graph
rng = np.random.default_rng()

# with open(info, 'w') as f:
#     print(ts, file=f)
#     print('\n', file=f)
#     print('Analysed times', file=f)
#     print(times, file=f)
#     print('\n', file=f)


#%% filtering
# freq_filt = ac.ts.get_allele_frequencies(
#     ts,
#     sample_sets=[ac.ts.individuals_at(ts, time=0, population=admix_pop)]
# )
# ts = ts.delete_sites(np.where(freq_filt[0] < 0.05)[0])

#%%
def ts_reps(files: list):
    for f in files:
        yield tskit.load(f)

results = []
for ts in ts_reps(files):
    results.append(
        ac.ts.analyze_trees(
            ts,
            times,
            n_samples,
            admix_pop,
            refs,
            alpha_mask,
            rng,
        )
    )


# with open(info, 'a') as f:
#     print('Total variance:', file=f)
#     print(totvar, file=f)
#     print('G:', file=f)
#     print(G, file=f)
#     print('G non-corrected:', file=f)
#     print(G_nc, file=f)
#     print('A\':', file=f)
#     print(Ap, file=f)
#     print('Total variance adjusted:', file=f)
#     print(totvar_adj, file=f)
#     print('\n', file=f)


#%% transform results
totvar = []
G = []
G_nc = []
Ap = []
Q = []
covmat_nc = []
covmat = []
for r in results:
    (t, gnc, g, a) =  ac.stats_from_matrices(
        r['covmat'],
        r['admix_cov'],
        r['drift_err'],
    )
    totvar.append(t)
    G_nc.append(gnc)
    G.append(g)
    Ap.append(a)
    Q.append(r['Q'])
    covmat_nc.append(r['covmat'])
    covmat.append(r['covmat'] - r['admix_cov'] - r['drift_err'])

totvar = np.array(totvar)
G_nc = np.array(G_nc)
G = np.array(G)
Ap = np.array(Ap)
Q = np.stack(Q)
covmat_nc = np.stack(covmat_nc)
covmat = np.stack(covmat)
# # convert to CIs
totvar_CI = ac.get_ci(totvar)
G_nc_CI = ac.get_ci(G_nc)
G_CI = ac.get_ci(G)
Ap_CI = ac.get_ci(Ap)
covmat_nc_CI = ac.get_ci(covmat_nc)
covmat_CI = ac.get_ci(covmat)

Q_CIs = [
    ac.get_ci(Q[:,:,i])
    for i in range(Q.shape[-1])
]


# %%
time_padding = 10

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

for i in range(Q.shape[2]):
    ac.plot_ci_line(x=times, CI=Q_CIs[i], ax=axs[0,0], color=colors_oi[i], label=f"Pop{i}")
axs[0,0].set_xlim(times[0] + time_padding, -time_padding)
axs[0,0].set_ylim((0,1))
axs[0,0].set_ylabel("Mean ancestry")
axs[0,0].set_xlabel("Time point")
axs[0,0].legend(loc="upper right")

combined_ci = ac.combine_covmat_CIs(covmat_CI, covmat_nc_CI)
scale_max = (
    np.max(np.abs([np.nanmin(combined_ci[1] - np.diag(np.diag(combined_ci[1]))),
    np.nanmax(combined_ci[1] - np.diag(np.diag(combined_ci[1])))]))
)
ac.plot_covmat_ci(combined_ci, axs[0, 1], scale_max)
axs[0,1].set_title('covariance matrix (raw lower, corrected upper)')

ymin, ymax = (
    np.min(np.concatenate([np.tril(covmat_nc_CI[0], k=-1), np.tril(covmat_CI[0], k=-1)])),
    np.max(np.concatenate([np.tril(covmat_nc_CI[2], k=-1), np.tril(covmat_CI[2], k=-1)])),
)
ac.cov_lineplot(times, covmat_nc_CI, axs[1, 0], colors=colors_oi, marker='o')
axs[1, 0].set_ylim(ymin, ymax)
axs[1, 0].set_ylabel('raw covariance (without bias)')
ac.cov_lineplot(times, covmat_CI, axs[1, 1], colors=colors_oi, marker='o')
axs[1, 1].set_ylim(ymin, ymax)
axs[1, 1].set_ylabel('admixture corrected covariance')

ac.plot_ci_line(x=times[1:], CI=totvar_CI, ax=axs[2, 0], marker='o')
axs[2, 0].set_xlim(times[1] + time_padding, -time_padding)
axs[2, 0].set_ylim(0)
axs[2, 0].set_ylabel('Total variance (t)')

ac.plot_ci_line(times[1:], G_nc_CI, ax=axs[2, 1], marker='o', linestyle='dashed')
ac.plot_ci_line(times[1:], G_CI, ax=axs[2, 1], marker='o')
ac.plot_ci_line(times[1:], Ap_CI, ax=axs[2, 1], marker='s', color='blue')
ymin, ymax = (
    np.min(np.concatenate([G_CI[0], G_nc_CI[0], Ap_CI[0]])),
    np.max(np.concatenate([G_CI[2], G_nc_CI[2], Ap_CI[2]])),
)
axs[2, 1].set_xlim(times[1] + time_padding, -time_padding)
axs[2, 1].set_ylim(ymin, ymax)
axs[2, 1].hlines(y=0, xmin=0, xmax=times[1] + time_padding, linestyles='dotted', colors='black')
axs[2, 1].set_xlabel('time')
axs[2, 1].set_ylabel("G(t) or A'(t)")

fig.tight_layout()

fig.savefig(
    snakemake.output['main_fig'],
)

#%%
if 'slim' in files[0]:
    fig, ax = plt.subplots(figsize=(6, 4))
    ztb = pd.read_csv(files[0].replace('.trees', '_pheno.tsv'), sep='\t')
    for f in files[1:]:
        ztb = pd.concat([ztb, pd.read_csv(f.replace('.trees', '_pheno.tsv'), sep='\t')])
    ztb['bgen'] = ztb.gen.max() - ztb.gen
    sns.lineplot(
        ztb[ztb.bgen < times[0] + 20], 
        x='bgen', y='mean_z', style='pop', hue='pop',
        estimator='mean', errorbar='ci', # 95% ci by default
        ax=ax,
    )
    ax.set_xlim(xmin=times[0] + 20, xmax=-5)
    ax.set_xlabel('generations')
    ax.set_ylabel('mean phenotype')
    fig.savefig(snakemake.output['main_fig'].replace('.pdf', '_pheno.pdf'))