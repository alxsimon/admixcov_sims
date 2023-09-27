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
        G_nde_CI,
        V_CI,
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

fig, axs = plt.subplots(3, 2, figsize=(10, 8), layout="tight")

# sci notation in colorbar
import matplotlib.ticker as tkr
formatter = tkr.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((0, 0))

k, l = (0, 0)
for i in range(len(Q_CIs)):
    ac.plot_ci_line(x=times, CI=Q_CIs[i], ax=axs[0,0], color=colors_oi[i], label=f"Pop{i}", marker='o')
axs[k, l].set_xlim(times[0] + time_padding, times[-1] - time_padding)
axs[k, l].set_ylim((0,1))
axs[k, l].set_ylabel("Mean ancestry proportion")
axs[k, l].set_xlabel("Time point")
axs[k, l].legend(loc="upper right")
axs[k, l].set_title('$A$', loc='left', fontdict={'fontweight': 'bold'})

k, l = (0, 1)
combined_ci = ac.combine_covmat_CIs(covmat_CI, covmat_nc_CI)
scale_max = (
    np.max(np.abs([np.nanmin(combined_ci[1] - np.diag(np.diag(combined_ci[1]))),
    np.nanmax(combined_ci[1] - np.diag(np.diag(combined_ci[1])))]))
)
ac.plot_covmat_ci(
    combined_ci,
    axs[k, l],
    scale_max,
    cbar_kws={'label': 'covariance', "format": formatter},
)
axs[k, l].set_title("B", loc='left', fontdict={'fontweight': 'bold'})

x_shift = 2
k, l = (1, 0)
ac.cov_lineplot(times, covmat_nc_CI, axs[k, l], colors=colors_oi, d=2)
axs[k, l].set_ylabel('Before correction')
axs[k, l].set_title("C", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].set_xlim(times[1] + x_shift, times[-2] - 4 * x_shift)
axs[k, l].hlines(y=0, xmin=times[1] + x_shift, xmax=times[-2] - 4 * x_shift, linestyles='dotted', colors='grey')
axs[k, l].yaxis.set_major_formatter(formatter)
k, l = (1, 1)
ac.cov_lineplot(times, covmat_CI, axs[k, l], colors=colors_oi, d=2, ylim=axs[1, 0].get_ylim())
axs[k, l].set_ylabel('After correction')
axs[k, l].set_title("D", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].set_xlim(times[1] + x_shift, times[-2] - 4 * x_shift)
axs[k, l].hlines(y=0, xmin=times[1] + x_shift, xmax=times[-2] - 4 * x_shift, linestyles='dotted', colors='grey')
axs[k, l].yaxis.set_major_formatter(formatter)

k, l = (2, 0)
ac.plot_ci_line(x=times[1:], CI=totvar_CI, ax=axs[k, l], marker='o')
axs[k, l].set_xlim(times[1] + time_padding, times[-1] - time_padding)
axs[k, l].set_ylim(0)
axs[k, l].set_ylabel('Total variance (t)')
axs[k, l].set_title("E", loc='left', fontdict={'fontweight': 'bold'})

k, l = (2, 1)
x_shift = 2
ymin = np.min([G_CI[1], G_nc_CI[1], G_nde_CI[1]]) - 0.1
ac.plot_ci_line(times[1:] + x_shift, G_nc_CI, ax=axs[k, l], marker='o', linestyle='dashed', label='$G_{nc}$')
ac.plot_ci_line(times[1:] + 2 * x_shift, G_nde_CI, ax=axs[k, l], marker='^', linestyle='dashdot', label='$G_{nde}$')
ac.plot_ci_line(times[1:], G_CI, ax=axs[k, l], marker='o', label='$G$')
ac.plot_ci_line(times[1:] - x_shift, Ap_CI, ax=axs[k, l], marker='s', color='blue', label='$A$')
axs[k, l].set_xlim(times[1] + time_padding, times[-1] - time_padding)
axs[k, l].set_ylim(ymax=1.1, ymin=ymin)
axs[k, l].hlines(y=0, xmin=times[-1] - time_padding, xmax=times[1] + time_padding, linestyles='dotted', colors='grey')
axs[k, l].set_xlabel('t')
axs[k, l].set_ylabel("Proportion of variance ($p_t - p_{160}$)")
axs[k, l].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[k, l].set_title("F", loc='left', fontdict={'fontweight': 'bold'})
for i, t in enumerate(times[1:]):
    (_, ytop) = axs[k, l].get_ylim()
    if G_CI[0][i]*G_CI[2][i] > 0:
        axs[k, l].annotate("*", xy=(t, ytop))

fig.savefig(
    snakemake.output['main_fig'],
)

#%%
if 'slim' in snakemake.input['pickle']:
    n_traits = 3
    fig, axs = plt.subplots(1, n_traits, figsize=(n_traits*4, 3), layout='constrained')
    for i in range(n_traits):
        sns.lineplot(
            ztb,
            x='bgen', y=f'mean_z{i}', style='pop', hue='pop',
            estimator='mean', errorbar='ci', # 95% ci by default
            ax=axs[i],
        )
        axs[i].set_xlim(xmin=200, xmax=-5)
        axs[i].set_xlabel('generations')
        axs[i].set_ylabel(f'mean phenotype {i}')
    fig.savefig(snakemake.output['pheno_fig'])