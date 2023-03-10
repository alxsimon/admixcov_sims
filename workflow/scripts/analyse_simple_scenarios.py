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
unit_n_sample = snakemake.params['n_sample']
unit_ref_n_sample = snakemake.params['ref_n_sample']
census_time = snakemake.params['census_time']
# info = snakemake.output['info']

drop_times = 2 if 'slim' in files[0] else 1

ts = tskit.load(files[0]) # extract info common to all trees
times = np.flip(ac.ts.get_times(ts))[drop_times:]
n_samples = [unit_n_sample]*len(times)
graph = demes.load(demes_file)
admix_pop = len(graph.demes) - 1
# not using ts.num_populations here as pyslim adds an additional one in ts
refs = [
    {'pop': i, 'time': census_time, 'n': unit_ref_n_sample}
    for i in range(admix_pop)
]
alpha_mask = np.array(
    [p.proportions for p in graph.pulses]
) > 0 # create alphas from graph
rng = np.random.default_rng()

# with open(info, 'w') as f:
#     print(ts, file=f)
#     print('\n', file=f)
#     print('Analysed times', file=f)
#     print(times, file=f)
#     print('\n', file=f)


#%% filtering
# freq_filt = ac.ts.get_allele_frequencies(
#     ts,
#     sample_sets=[ac.ts.individuals_at(ts, time=0, population=admix_pop)]
# )
# ts = ts.delete_sites(np.where(freq_filt[0] < 0.05)[0])

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
            admix_pop,
            refs,
            alpha_mask,
            rng,
        )
    )


# with open(info, 'a') as f:
#     print('Total variance:', file=f)
#     print(totvar, file=f)
#     print('G:', file=f)
#     print(G, file=f)
#     print('G non-corrected:', file=f)
#     print(G_nc, file=f)
#     print('A\':', file=f)
#     print(Ap, file=f)
#     print('Total variance adjusted:', file=f)
#     print(totvar_adj, file=f)
#     print('\n', file=f)


#%% transform results
totvar = []
G = []
G_nc = []
Ap = []
Q = []
covmat_nc = []
covmat = []
for r in results:
    (t, gnc, g, a) =  ac.stats_from_matrices(
        r['covmat'],
        r['admix_cov'],
        r['drift_err'],
    )
    totvar.append(t)
    G_nc.append(gnc)
    G.append(g)
    Ap.append(a)
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
# # convert to CIs
totvar_CI = ac.get_ci(totvar)
G_nc_CI = ac.get_ci(G_nc)
G_CI = ac.get_ci(G)
Ap_CI = ac.get_ci(Ap)
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


with open(snakemake.output['pickle'], 'wb') as fw:
    pickle.dump(
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
        ),
        fw
    )