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

prefix_gss = snakemake.params['prefix_gss']
prefix_bgs = snakemake.params['prefix_bgs']
recombinations = snakemake.params['rec']
intervals = snakemake.params['intervals']
res = {k: {} for k in recombinations}
sum_cov = {}
sum_var = {}
sum_varcov = {}

fig, axs = plt.subplots(1, 3, figsize=(12, 3.5), layout="constrained")

# GSS plots on axes 0 and 2
for rec in recombinations:
    for i in intervals:
        with open(f'{prefix_gss}_r{rec}_i{i}.pickle', 'rb') as fr:
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
    ac.plot_ci_line(x_val, sum_var[rec], axs[2], color=colors_oi[j], label=f"r={rec}", marker='o')
    ymin = min(ymin, np.min(np.stack([sum_cov[rec][0], sum_var[rec][0]])))
    ymax = max(ymax, np.max(np.stack([sum_cov[rec][2], sum_var[rec][2]])))

ylim = (ymin, ymax)
axs[0].legend(ncols=2)
axs[0].set_xlabel("sampling interval (gen)")
axs[0].set_ylabel("Sum(Cov)")
_ = axs[0].set_ylim(ylim)
axs[2].set_xlabel("sampling interval (gen)")
axs[2].set_ylabel("Sum(Var)")
_ = axs[2].set_ylim(ylim)

_ = axs[0].set_xticks([int(x) + 0.5 for x in intervals])
_ = axs[0].set_xticklabels(intervals)
_ = axs[2].set_xticks([int(x) + 0.5 for x in intervals])
_ = axs[2].set_xticklabels(intervals)

# a supplementary figure for just r=2e-8
supfig, supaxs = plt.subplots(1, 2, figsize=(10, 3.5), layout="constrained")
x_val = [int(x) for x in intervals]
rec = "2e-8"
ac.plot_ci_line(x_val, sum_cov[rec], supaxs[0], marker='o', color=colors_oi[3])
ac.plot_ci_line(x_val, sum_var[rec], supaxs[1], marker='o', color=colors_oi[3])
ylim = (
    np.min(np.stack([sum_cov[rec][0], sum_var[rec][0]])),
    np.max(np.stack([sum_cov[rec][2], sum_var[rec][2]]))
)
supaxs[0].set_xlabel("sampling interval (gen)")
supaxs[0].set_ylabel("Sum(Cov)")
supaxs[0].set_title("A", loc='left', fontdict={'fontweight': 'bold'})
supaxs[0].set_title('Covariances GSS')
_ = supaxs[0].set_xticks(x_val)
_ = supaxs[0].set_xticklabels(intervals)
_ = supaxs[0].set_ylim(ylim)
supaxs[1].set_xlabel("sampling interval (gen)")
supaxs[1].set_ylabel("Sum(Var)")
supaxs[1].set_title("B", loc='left', fontdict={'fontweight': 'bold'})
supaxs[1].set_title('Variances GSS')
_ = supaxs[1].set_xticks(x_val)
_ = supaxs[1].set_xticklabels(intervals)
_ = supaxs[1].set_ylim(ylim)
supaxs[0].yaxis.set_major_formatter(formatter)
supaxs[1].yaxis.set_major_formatter(formatter)
supfig.savefig(snakemake.output['supfig'])

# BGS plot on axes 1
for rec in recombinations:
    for i in intervals:
        with open(f'{prefix_bgs}_r{rec}_i{i}.pickle', 'rb') as fr:
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

# ymin = np.Inf
# ymax = -np.Inf
for j, rec in enumerate(recombinations):
    x_val = [int(x) + 0.25 * j for x in intervals]
    ac.plot_ci_line(x_val, sum_cov[rec], axs[1], label=f"r={rec}", marker='o', color=colors_oi[j])
#     ymin = min(ymin, np.min(sum_cov[rec][0]))
#     ymax = max(ymax, np.max(sum_cov[rec][2]))

ylim = (ymin, ymax)
axs[1].set_xlabel("sampling interval (gen)")
axs[1].set_ylabel("Sum(Cov)")
_ = axs[1].set_ylim(ymin=0)

_ = axs[1].set_xticks([int(x) + 0.5 for x in intervals])
_ = axs[1].set_xticklabels(intervals)

# ==============
axs[0].set_title("A", loc='left', fontdict={'fontweight': 'bold'})
axs[0].set_title('Covariances GSS')

axs[1].set_title("B", loc='left', fontdict={'fontweight': 'bold'})
axs[1].set_title('Covariances BGS')

axs[2].set_title("C", loc='left', fontdict={'fontweight': 'bold'})
axs[2].set_title('Variances GSS')

axs[0].yaxis.set_major_formatter(formatter)
axs[1].yaxis.set_major_formatter(formatter)
axs[2].yaxis.set_major_formatter(formatter)

fig.savefig(snakemake.output['fig'])