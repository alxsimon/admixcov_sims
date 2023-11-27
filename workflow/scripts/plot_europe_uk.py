#%%
import admixcov as ac
import demes
import demesdraw
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

demes_file = snakemake.input['demes_file']
graph = demes.load(demes_file)
fig, ax = plt.subplots(figsize=(8, 8))
demesdraw.tubes(graph, log_time=True, ax=ax)
fig.savefig(snakemake.output['fig_demo'])

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

import matplotlib.ticker as plticker
loc = plticker.MultipleLocator(base=1.0)

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
delta_list = [f"$\\Delta p_{{{int(t)}}}$" for t in range(len(times) - 1)]

# sci notation in colorbar
import matplotlib.ticker as tkr
formatter = tkr.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((0, 0))

fig, axs = plt.subplots(2, 2, figsize=(10, 8), layout='constrained')

k, l = (0, 1)
fmts = ['-o', '-s', '-^']
labels = ['WHG', 'ANA', 'YAM']
for i, pop in enumerate(labels):
    ac.plot_ci_line(x=times, CI=Q_CIs[i], ax=axs[k, l], color=colors_oi[i], label=pop, fmt=fmts[i])
for x1, x2, txt in zip(times[:-1], times[1:], delta_list):
    _ = axs[k, l].text(x2+(x1 - x2)/2, 0.9, txt, ha='center')
for i, t in enumerate(times):
    _ = axs[k, l].text(t, 0.8, str(i), ha='center')
for x1, x2 in zip(times[1::2], times[2::2]):
    _ = axs[k, l].axvspan(x1, x2, facecolor='grey', alpha=0.10)
axs[k, l].set_xlim(times[0] + time_padding, times[-1] - time_padding)
axs[k, l].set_ylim(top=1)
axs[k, l].set_ylabel("Mean ancestry proportion")
axs[k, l].set_xlabel("Time (years BP)")
axs[k, l].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
axs[k, l].set_title("B", loc='left', fontdict={'fontweight': 'bold'})

x_shift = 0.1
new_times = np.array(range(len(times)))
k, l = (0, 0)
ac.cov_lineplot(new_times, covmat_nc_CI, axs[k, l], colors=colors_oi, d=x_shift, labels=delta_list)
axs[k, l].set_xlim(new_times[1] - x_shift, new_times[-2] + 3 * x_shift)
axs[k, l].hlines(y=0, xmin=0, xmax=new_times[-1] + 3 * x_shift, linestyles='dotted', colors='grey')
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_t$)")
axs[k, l].set_xlabel("t")
axs[k, l].set_title('Before admix. correction')
axs[k, l].set_title("A", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].xaxis.set_major_locator(loc)
axs[k, l].yaxis.set_major_formatter(formatter)

k, l = (1, 0)
ac.cov_lineplot(new_times, covmat_CI, axs[k, l], colors=colors_oi, d=x_shift, labels=delta_list, ylim=axs[0, 0].get_ylim())
axs[k, l].set_xlim(new_times[1] - x_shift, new_times[-2] + 3 * x_shift)
axs[k, l].hlines(y=0, xmin=0, xmax=new_times[-1] + 3 * x_shift, linestyles='dotted', colors='grey')
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_t$)")
axs[k, l].set_xlabel('t')
axs[k, l].set_title('After admix. correction')
axs[k, l].set_title("C", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), title="$\\Delta p_i$", ncol=3)
axs[k, l].xaxis.set_major_locator(loc)
axs[k, l].yaxis.set_major_formatter(formatter)

k, l = (1, 1)
ac.plot_ci_line(new_times[1:] + x_shift, G_nc_CI, ax=axs[k, l], linestyle='dashed', marker='o', label='$G_{nc}$')
ac.plot_ci_line(new_times[1:], G_CI, ax=axs[k, l], marker='o', label='$G$')
ac.plot_ci_line(new_times[1:] - x_shift, Ap_CI, ax=axs[k, l], color='blue', marker='s', label='$A$')
axs[k, l].set_xlim(new_times[1] - 2*x_shift, new_times[-1] + 2*x_shift)
axs[k, l].hlines(y=0, xmin=new_times[-1], xmax=new_times[1], colors='grey', linestyles='dotted')
axs[k, l].set_ylim(ymax=1)
axs[k, l].set_xlabel('t')
axs[k, l].set_ylabel("Proportion of variance ($p_t - p_{0}$)")
axs[k, l].legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
axs[k, l].set_title("D", loc='left', fontdict={'fontweight': 'bold'})
axs[k, l].xaxis.set_major_locator(loc)
for i, t in enumerate(new_times[1:]):
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