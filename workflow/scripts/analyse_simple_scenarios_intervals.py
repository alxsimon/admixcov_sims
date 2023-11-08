#%%
import admixcov as ac
import tskit
import demes
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle

#%%
files = snakemake.input['files']
demes_file = snakemake.input['demes_file']
census_time = snakemake.params['census_time']
interval = int(snakemake.wildcards['inter'])
start_sampling = snakemake.params['start_sampling']

# drop_times = 2 if 'slim' in files[0] else 1

times = np.flip(range(0, start_sampling + 1, interval))

ts = tskit.load(files[0]) # extract info common to all trees
# times = np.flip(ac.ts.get_times(ts))[drop_times:]
graph = demes.load(demes_file)
N_admix_pop = len(graph.demes) - 1

n_samples = [snakemake.params['n_sample']]*len(times)
assert len(n_samples) == len(times)

ref_n_sample = snakemake.params['ref_n_sample']

# not using ts.num_populations here as pyslim adds an additional one in ts
refs = [
    {'pop': i, 'time': census_time, 'n': ref_n_sample}
    for i in range(N_admix_pop)
]

proportions = np.zeros((times.size - 1, len(graph.pulses[0].proportions)))
alpha_times = np.array([p.time for p in graph.pulses])
graph_proportions = np.array(
    [p.proportions for p in graph.pulses]
)
for i, t in enumerate(times[:-1]):
    where = np.where((alpha_times >= times[i+1]) & (alpha_times < times[i]))[0]
    if where.size > 0:
        proportions[i] = graph_proportions[where]
# This finds where to put the pulses when we have more time steps than pulses.
# We take into account that slim sampling is done after migration, so that
# a pulse will be taken into account into the previous time interval
# and not the following.

alpha_mask = np.array(
    proportions
) > 0 # create alphas from graph
rng = np.random.default_rng()

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
            N_admix_pop,
            refs,
            alpha_mask,
            rng,
        )
    )

#%% transform results
totvar = []
G = []
G_nc = []
G_nde = []
Ap = []
V = []
Q = []
covmat_nc = []
covmat = []
for r in results:
    (t, gnc, g, a, gnde, v) =  ac.stats_from_matrices(
        r['covmat'],
        r['admix_cov'],
        r['drift_err'],
    )
    totvar.append(np.array(t) / r['hz'][0]) # dividing by first time point hz
    G_nc.append(gnc)
    G.append(g)
    Ap.append(a)
    G_nde.append(gnde)
    V.append(v)
    Q.append(r['Q'])
    covmat_nc.append(r['covmat'])
    covmat.append(r['covmat'] - r['admix_cov'] - r['drift_err'])

totvar = np.array(totvar)
G_nc = np.array(G_nc)
G = np.array(G)
Ap = np.array(Ap)
G_nde = np.array(G_nde)
V = np.array(V)
Q = np.stack(Q)
covmat_nc = np.stack(covmat_nc)
covmat = np.stack(covmat)

sum_varcov = covmat.sum(axis=(1, 2))
sum_var = covmat.diagonal(0, 1, 2).sum(axis=1)
sum_cov = sum_varcov - sum_var

# convert to CIs
totvar_CI = ac.get_ci(totvar)
G_nc_CI = ac.get_ci(G_nc)
G_CI = ac.get_ci(G)
Ap_CI = ac.get_ci(Ap)
G_nde_CI = ac.get_ci(G_nde)
V_CI = ac.get_ci(V)

covmat_nc_CI = ac.get_ci(covmat_nc)
covmat_CI = ac.get_ci(covmat)

sum_varcov_CI = ac.get_ci(sum_varcov)
sum_var_CI = ac.get_ci(sum_var)
sum_cov_CI = ac.get_ci(sum_cov)

Q_CIs = [
    ac.get_ci(Q[:,:,i])
    for i in range(Q.shape[-1])
]

if 'slim' in files[0]:
    ztb = pd.read_csv(files[0].replace('.trees', '_pheno.tsv'), sep='\t')
    for f in files[1:]:
        ztb = pd.concat([ztb, pd.read_csv(f.replace('.trees', '_pheno.tsv'), sep='\t')])
    ztb['bgen'] = ztb.gen.max() - ztb.gen
else:
    ztb = None

with open(snakemake.output['pickle'], 'wb') as fw:
    pickle.dump(
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
            sum_varcov_CI,
            sum_var_CI,
            sum_cov_CI,
            Q_CIs,
            ztb,
        ),
        fw
    )