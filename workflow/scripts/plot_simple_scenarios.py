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
        G_de_CI,
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

times = np.array(times) # ensure it is an array

fig, axs = plt.subplots(3, 2, figsize=(10, 8))

k, l = (0, 0)
for i in range(len(Q_CIs)):
    ac.plot_ci_line(x=times, CI=Q_CIs[i], ax=axs[0,0], color=colors_oi[i], label=f"Pop{i}", marker='o')
axs[k, l].set_xlim(times[0] + time_padding, times[-1] - time_padding)
# axs[k, l].set_ylim((0,1))
axs[k, l].set_ylabel("Mean ancestry")
axs[k, l].set_xlabel("Time point")
axs[k, l].legend(loc="upper right")

k, l = (0, 1)
combined_ci = ac.combine_covmat_CIs(covmat_CI, covmat_nc_CI)
scale_max = (
    np.max(np.abs([np.nanmin(combined_ci[1] - np.diag(np.diag(combined_ci[1]))),
    np.nanmax(combined_ci[1] - np.diag(np.diag(combined_ci[1])))]))
)
ac.plot_covmat_ci(combined_ci, axs[k, l], scale_max)
axs[k, l].set_title('covariance matrix (raw lower, corrected upper)')

k, l = (1, 0)
ac.cov_lineplot(times, covmat_nc_CI, axs[k, l], colors=colors_oi, marker='o', time_padding=time_padding, d=2)
axs[k, l].set_ylabel('raw covariance (without bias)')
k, l = (1, 1)
ac.cov_lineplot(times, covmat_CI, axs[k, l], colors=colors_oi, marker='o', time_padding=time_padding, d=2, ylim=axs[k, l].get_ylim())
axs[k, l].set_ylabel('admixture corrected covariance')

k, l = (2, 0)
ac.plot_ci_line(x=times[1:], CI=totvar_CI, ax=axs[k, l], marker='o')
axs[k, l].set_xlim(times[1] + time_padding, times[-1] - time_padding)
axs[k, l].set_ylim(0)
axs[k, l].set_ylabel('Total variance (t)')

k, l = (2, 1)
x_shift = 2
ac.plot_ci_line(times[1:] + x_shift, G_nc_CI, ax=axs[k, l], marker='o', linestyle='dashed')
ac.plot_ci_line(times[1:] + 2 * x_shift, G_de_CI, ax=axs[k, l], marker='^', linestyle='dotted')
ac.plot_ci_line(times[1:], G_CI, ax=axs[k, l], marker='o')
ac.plot_ci_line(times[1:] - x_shift, Ap_CI, ax=axs[k, l], marker='s', color='blue')
axs[k, l].set_xlim(times[1] + time_padding, times[-1] - time_padding)
axs[k, l].hlines(y=0, xmin=times[-1] - time_padding, xmax=times[1] + time_padding, linestyles='dotted', colors='black')
axs[k, l].set_xlabel('time')
axs[k, l].set_ylabel("G(t) or A'(t)")
for i, t in enumerate(times[1:]):
    if G_CI[0][i]*G_CI[2][i] > 0:
        axs[k, l].annotate("*", xy=(t, 0.1))

fig.tight_layout()

fig.savefig(
    snakemake.output['main_fig'],
)

#%%
if 'slim' in snakemake.input['pickle']:
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.lineplot(
        # ztb[ztb.bgen < times[0] + 20], 
        ztb,
        x='bgen', y='mean_z', style='pop', hue='pop',
        estimator='mean', errorbar='ci', # 95% ci by default
        ax=ax,
    )
    ax.set_xlim(xmin=200, xmax=-5)
    ax.set_xlabel('generations')
    ax.set_ylabel('mean phenotype')
    fig.savefig(snakemake.output['pheno_fig'])