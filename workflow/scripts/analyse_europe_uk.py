#%%
import admixcov as ac
import tskit
import demes
import demesdraw
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pickle

#%%
files = snakemake.input['files']
demes_file = snakemake.input['demes_file']
census_time = snakemake.params['census_time']

# drop_times = 2 if 'slim' in files[0] else 1

ts = tskit.load(files[0]) # extract info common to all trees
# times = np.flip(ac.ts.get_times(ts))[drop_times:]
times = [150, 130, 110, 90, 70, 50, 0]
graph = demes.load(demes_file)

fig, ax = plt.subplots(figsize=(8, 8))
demesdraw.tubes(graph, log_time=True, ax=ax)
fig.savefig(snakemake.output['fig_demo'])

n_samples = snakemake.params['n_samples']
assert len(n_samples) == len(times)

ref_n_samples = snakemake.params['ref_n_samples']

# WHG, ANA, YAM
refs = [
    {'pop': p, 'time': c, 'n': n}
    for (p, n, c) in zip([5, 4, 7], ref_n_samples, [200, 200, 150])
]
alpha_mask = np.array([ # WHG, ANA, YAM
    [0, 0, 1],
    [0, 1, 0],
    [0, 1, 0],
    [0, 1, 0],
    [1, 0, 0],
    [0, 1, 0],
], dtype=bool)
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
            8, # focal pop
            refs,
            alpha_mask,
            rng,
        )
    )

#%% transform results
totvar = []
G = []
G_nc = []
Ap = []
G_de = []
Q = []
covmat_nc = []
covmat = []
for r in results:
    (t, gnc, g, a, gde) =  ac.stats_from_matrices(
        r['covmat'],
        r['admix_cov'],
        r['drift_err'],
    )
    totvar.append(np.array(t) / r['hz'][0]) # dividing by first time point hz
    G_nc.append(gnc)
    G.append(g)
    Ap.append(a)
    G_de.append(gde)
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
# convert to CIs
totvar_CI = ac.get_ci(totvar)
G_nc_CI = ac.get_ci(G_nc)
G_CI = ac.get_ci(G)
Ap_CI = ac.get_ci(Ap)
G_de_CI = ac.get_ci(G_de)

covmat_nc_CI = ac.get_ci(covmat_nc)
covmat_CI = ac.get_ci(covmat)

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
            G_de_CI,
            covmat_nc_CI,
            covmat_CI,
            Q_CIs,
            ztb,
        ),
        fw
    )