import admixcov as ac
import numpy as np
import matplotlib.pyplot as plt
import pickle
# sci notation formatter
import matplotlib.ticker as tkr
formatter = tkr.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((0, 0))

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

# Plot selection scenario
prefix = snakemake.params['prefix']
recombinations = snakemake.params['rec']
intervals = snakemake.params['intervals']
res = {k: {} for k in recombinations}
sum_cov = {}
sum_var = {}
sum_varcov = {}

fig, axs = plt.subplots(1, 2, figsize=(12, 4), layout="constrained")

for rec in recombinations:
    for i in intervals:
        with open(f'{prefix}_r{rec}_i{i}.pickle', 'rb') as fr:
            res[rec][i] = {}
            (
                res[rec][i]['times'],
                res[rec][i]['totvar_CI'],
                res[rec][i]['G_nc_CI'],
                res[rec][i]['G_CI'],
                res[rec][i]['Ap_CI'],
                res[rec][i]['G_nde_CI'],
                res[rec][i]['V_CI'],
                res[rec][i]['covmat_nc_CI'],
                res[rec][i]['covmat_CI'],
                res[rec][i]['sum_varcov_CI'],
                res[rec][i]['sum_var_CI'],
                res[rec][i]['sum_cov_CI'],
                res[rec][i]['Q_CIs'],
                res[rec][i]['ztb'],
            ) = pickle.load(fr)

    sum_cov[rec] = tuple([
        np.array([res[rec][x]['sum_cov_CI'][i] for x in intervals])
        for i in range(3)
    ])
    sum_var[rec] = tuple([
        np.array([res[rec][x]['sum_var_CI'][i] for x in intervals])
        for i in range(3)
    ])

ymin = np.Inf
ymax = -np.Inf
for j, rec in enumerate(recombinations):
    x_val = [int(x) + 0.25 * j for x in intervals]
    ac.plot_ci_line(x_val, sum_cov[rec], axs[0], label=f"r={rec}", marker='o', color=colors_oi[j])
    ac.plot_ci_line(x_val, sum_var[rec], axs[1], color=colors_oi[j], label=f"r={rec}", marker='o')
    ymin = min(ymin, np.min(np.stack([sum_cov[rec][0], sum_var[rec][0]])))
    ymax = max(ymax, np.max(np.stack([sum_cov[rec][2], sum_var[rec][2]])))

ylim = (ymin, ymax)
axs[0].legend(ncols=3)
axs[0].set_xlabel("sampling interval (gen)")
axs[0].set_ylabel("Sum(Cov)")
axs[0].set_title("A", loc='left', fontdict={'fontweight': 'bold'})
axs[0].set_title(f"Covariances {snakemake.params['model']}")
_ = axs[0].set_ylim(ylim)
axs[1].set_xlabel("sampling interval (gen)")
axs[1].set_ylabel("Sum(Var)")
axs[1].set_title("B", loc='left', fontdict={'fontweight': 'bold'})
axs[1].set_title(f"Variances {snakemake.params['model']}")
_ = axs[1].set_ylim(ylim)

_ = axs[0].set_xticks([int(x) + 0.5 for x in intervals])
_ = axs[0].set_xticklabels(intervals)
_ = axs[1].set_xticks([int(x) + 0.5 for x in intervals])
_ = axs[1].set_xticklabels(intervals)
axs[0].yaxis.set_major_formatter(formatter)
axs[1].yaxis.set_major_formatter(formatter)

fig.savefig(snakemake.output['fig'])