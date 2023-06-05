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

fig, axs = plt.subplots(2, 2, figsize=(10, 8), layout='constrained')

k, l = (0, 1)
fmts = ['-o', '-s', '-^']
labels = ['WHG', 'ANA', 'YAM']
for i, pop in enumerate(labels):
    ac.plot_ci_line(x=times, CI=Q_CIs[i], ax=axs[k, l], color=colors_oi[i], label=pop, fmt=fmts[i])
axs[k, l].set_xlim(times[0] + time_padding, times[-1] - time_padding)
axs[k, l].set_ylim(top=1)
axs[k, l].set_ylabel("Mean ancestry")
axs[k, l].set_xlabel("Time (generations BP)")
axs[k, l].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
axs[k, l].set_title("B", loc='left', fontdict={'fontweight': 'bold'})

# combined_ci = ac.combine_covmat_CIs(covmat_CI, covmat_nc_CI)
# scale_max = (
#     np.max(np.abs([np.nanmin(combined_ci[1] - np.diag(np.diag(combined_ci[1]))),
#     np.nanmax(combined_ci[1] - np.diag(np.diag(combined_ci[1])))]))
# )
# ac.plot_covmat_ci(combined_ci, axs[0, 1], scale_max)
# axs[0,1].set_title('covariance matrix (raw lower, corrected upper)')

x_shift = 2
k, l = (0, 0)
ac.cov_lineplot(times, covmat_nc_CI, axs[k, l], colors=colors_oi, time_padding=time_padding, d=x_shift, marker='o')
axs[k, l].set_xlim(times[1] + x_shift + time_padding, times[-2] - x_shift - time_padding)
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_t$)")
axs[k, l].set_xlabel("t")
axs[k, l].set_title('Before admix. correction')
axs[k, l].set_title("A", loc='left', fontdict={'fontweight': 'bold'})

k, l = (1, 0)
ac.cov_lineplot(times, covmat_CI, axs[k, l], colors=colors_oi, time_padding=time_padding, d=x_shift, marker='o', ylim=axs[0, 0].get_ylim())
axs[k, l].set_xlim(times[1] + x_shift + time_padding, times[-2] - x_shift - time_padding)
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_t$)")
axs[k, l].set_xlabel('t')
axs[k, l].set_title('After admix. correction')
axs[k, l].set_title("C", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), title="$\\Delta p_i$", ncol=3)

# ac.plot_ci_line(x=times[1:], CI=totvar_CI, ax=axs[2, 0], marker='o')
# axs[2, 0].set_xlim(times[1] + time_padding, times[-1] - time_padding)
# axs[2, 0].set_ylim(0)
# axs[2, 0].set_ylabel('Total variance (t)')

k, l = (1, 1)
ac.plot_ci_line(times[1:] + x_shift, G_nc_CI, ax=axs[k, l], linestyle='dashed', marker='o', label='$G_{nc}$')
ac.plot_ci_line(times[1:] + 2 * x_shift, G_de_CI, ax=axs[k, l], marker='^', linestyle='dashdot', label='$G_{de}$')
ac.plot_ci_line(times[1:], G_CI, ax=axs[k, l], marker='o', label='G')
ac.plot_ci_line(times[1:] - x_shift, Ap_CI, ax=axs[k, l], color='blue', marker='s', label='A\'')
axs[k, l].set_xlim(times[1] + x_shift + time_padding, times[-1] - x_shift - time_padding)
axs[k, l].hlines(y=0, xmin=times[-1] - time_padding, xmax=times[1] + time_padding, colors='grey', linestyles='dotted')
axs[k, l].set_xlabel('t')
axs[k, l].set_ylabel("Proportion of variance ($p_t - p_{5424}$)")
axs[k, l].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[k, l].set_title("D", loc='left', fontdict={'fontweight': 'bold'})
for i, t in enumerate(times[1:]):
    if G_CI[0][i]*G_CI[2][i] > 0:
        axs[k, l].annotate("*", xy=(t, 0.1))

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