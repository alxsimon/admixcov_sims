#%%
import admixcov as ac
import tskit
import demes
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle

with open(snakemake.input['pickle'], 'rb') as fr:
	(
		times,
        totvar_CI,
        G_nc_CI,
        G_CI,
        Ap_CI,
        covmat_nc_CI,
        covmat_CI,
        Q_CIs,
        ztb,	
	) = pickle.load(fr)

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
if 'slim' in snakemake.input['pickle']:
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.lineplot(
        ztb[ztb.bgen < times[0] + 20], 
        x='bgen', y='mean_z', style='pop', hue='pop',
        estimator='mean', errorbar='ci', # 95% ci by default
        ax=ax,
    )
    ax.set_xlim(xmin=200, xmax=-5)
    ax.set_xlabel('generations')
    ax.set_ylabel('mean phenotype')
    fig.savefig(snakemake.output['pheno_fig'])