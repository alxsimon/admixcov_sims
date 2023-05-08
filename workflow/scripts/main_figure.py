import admixcov as ac
import tskit
import demes
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle

fig, axs = plt.subplots(3, 2, figsize=(12, 10))

with open(snakemake.input['sim_neutral'], 'rb') as fr:
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
delta_list = [f"$\\Delta p_{{{int(t)}}}$" for t in times[:-1]]

k, l = (0, 0)
for i in range(len(Q_CIs)):
    ac.plot_ci_line(x=times, CI=Q_CIs[i], ax=axs[k, l], color=colors_oi[i], label=f"Pop{i}", marker='o')
for x, txt in zip([t - 10 for t in times[:-1]], delta_list):
	_ = axs[k, l].text(x, 1, txt, ha='center')
for x in times[1::2]:
    _ = axs[k, l].axvspan(x, x + 20, facecolor='grey', alpha=0.10)
for x in [110, 90, 30, 10]:
	_ = axs[k, l].annotate("", xy=(x, 0.1), xytext=(x, 0), arrowprops=dict(arrowstyle="->"))
axs[k, l].set_xlim(times[0] + time_padding, times[-1] - time_padding)
axs[k, l].set_ylabel("Mean ancestry")
axs[k, l].set_xlabel("Time (gen. BP)")
axs[k, l].legend(loc="center left")
axs[k, l].set_title("A", loc='left', fontdict={'fontweight': 'bold'})

k, l = (0, 1)
combined_ci = ac.combine_covmat_CIs(covmat_CI, covmat_nc_CI)
scale_max = (
    np.max(np.abs([np.nanmin(combined_ci[1] - np.diag(np.diag(combined_ci[1]))),
    np.nanmax(combined_ci[1] - np.diag(np.diag(combined_ci[1])))]))
)
ac.plot_covmat_ci(
	combined_ci, axs[k, l], scale_max,
	delta_labels=delta_list,
	cbar_kws={'label': 'covariance'},
)
axs[k, l].set_title("B", loc='left', fontdict={'fontweight': 'bold'})


k, l = (1, 0)
ac.cov_lineplot(times, covmat_nc_CI, axs[k, l], colors=colors_oi, marker='o', time_padding=time_padding, d=2)
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_t$)")
axs[k, l].set_xlabel('t')
axs[k, l].set_title('Before admix. correction')
axs[k, l].set_title("C", loc='left', fontdict={'fontweight': 'bold'})
k, l = (1, 1)
ac.cov_lineplot(times, covmat_CI, axs[k, l], colors=colors_oi, marker='o', time_padding=time_padding, d=2, ylim=axs[k, l - 1].get_ylim())
axs[k, l].set_ylabel("Cov($\\Delta p_i$, $\\Delta p_t$)")
axs[k, l].set_xlabel('t')
axs[k, l].set_title('After admix. correction')
axs[k, l].legend(loc='center left', bbox_to_anchor=(1, 0.5), title="$\\Delta p_i$")
axs[k, l].set_title("D", loc='left', fontdict={'fontweight': 'bold'})


k, l = (2, 0)
x_shift = 2
ac.plot_ci_line(times[1:] + 2, G_nc_CI, ax=axs[k, l], marker='o', linestyle='dashed', label='G_nc')
ac.plot_ci_line(times[1:], G_CI, ax=axs[k, l], marker='o', label='G')
ac.plot_ci_line(times[1:] - 2, Ap_CI, ax=axs[k, l], marker='s', color='blue', label='A\'')
axs[k, l].set_xlim(times[1] + time_padding, times[-1] - time_padding)
axs[k, l].hlines(y=0, xmin=times[-1] - time_padding, xmax=times[1] + time_padding, linestyles='dotted', colors='black')
axs[k, l].set_xlabel('t')
axs[k, l].set_ylabel("Proportion of variance ($p_t - p_{120}$)")
axs[k, l].set_title('neutral')
axs[k, l].set_title("E", loc='left', fontdict={'fontweight': 'bold'})

# ==================
with open(snakemake.input['sim_sel'], 'rb') as fr:
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


k, l = (2, 1)
x_shift = 2
ac.plot_ci_line(times[1:] + 2, G_nc_CI, ax=axs[k, l], marker='o', linestyle='dashed', label='G_nc')
ac.plot_ci_line(times[1:], G_CI, ax=axs[k, l], marker='o', label='G')
ac.plot_ci_line(times[1:] - 2, Ap_CI, ax=axs[k, l], marker='s', color='blue', label='A\'')
axs[k, l].set_xlim(times[1] + time_padding, times[-1] - time_padding)
axs[k, l].hlines(y=0, xmin=times[-1] - time_padding, xmax=times[1] + time_padding, linestyles='dotted', colors='black')
axs[k, l].set_xlabel('t')
axs[k, l].set_ylabel("Proportion of variance ($p_t - p_{120}$)")
axs[k, l].set_title('sel. jump t140')
axs[k, l].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs[k, l].set_title("F", loc='left', fontdict={'fontweight': 'bold'})

fig.tight_layout()
fig.savefig(snakemake.output['fig'])