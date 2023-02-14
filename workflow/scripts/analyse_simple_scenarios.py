import admixcov as ac
import tskit
import demes
import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
import pickle

trees_file = snakemake.input['trees_file']
demes_file = snakemake.input['demes_file']
unit_n_sample = snakemake.params['n_sample']
unit_ref_n_sample = snakemake.params['ref_n_sample']
census_time = snakemake.params['census_time']

ts = tskit.load(trees_file)
times = np.flip(ac.ts.get_times(ts))[1:]
n_samples = [unit_n_sample]*len(times)
graph = demes.load(demes_file)
admix_pop = ts.num_populations - 1
refs = [
    {'pop': i, 'time': census_time, 'n': unit_ref_n_sample}
    for i in range(admix_pop)
]
rng = np.random.default_rng()

samples_nodes = ac.ts.draw_sample_sets(ts, times, rng, admix_pop, n_samples)
ref_nodes = [
    ac.ts.draw_sample_set_at(ts, r['time'], rng, r['pop'], r['n'])
    for r in refs
]

ref_geno_sets, geno_sets = ac.ts.get_pseudohap_genotypes(
    ts, rng,
    samples=(ref_nodes, samples_nodes),
)

af = np.stack([
    np.mean(G, axis=1)
    for G in geno_sets
], axis=0)
ref_af = np.stack([
    np.mean(G, axis=1)
    for G in ref_geno_sets
], axis=0)

# admixture computations
samples_inds = [np.unique([ts.node(u).individual for u in s]) for s in samples_nodes]
ancestral_census_nodes = [
    np.where(
        (ts.tables.nodes.population == r['pop'])
        & (ts.tables.nodes.time == r['time'])
    )[0]
    for r in refs
]
Q = np.stack([
    np.mean(
        ac.ts.get_admixture_proportions(ts, g, ancestral_census_nodes),
        axis=0,
    )
    for g in samples_inds
])

alpha_mask = np.array(
    [p.proportions for p in graph.pulses]
) > 0 # create alphas from graph
alphas = ac.q2a_simple(Q, alpha_mask)
alphas[alphas < 0] = 0 # issues due to drift
covmat = ac.get_covariance_matrix(
    af, bias=True,
    sample_size=np.array(n_samples),
)
admix_cov = ac.get_admix_covariance_matrix(
    Q, ref_af, bias=True,
    ref_sample_size=np.array([d['n'] for d in refs]),
)
var_drift = ac.solve_for_variances(
    np.diag(covmat - admix_cov), alphas,
)
drift_err = ac.get_drift_err_matrix(
    var_drift, alphas,
)

G = ac.get_summary(
    covmat - admix_cov - drift_err,
    af, np.array(n_samples),
    include_diag=False, abs=False
)

Ap = ac.get_summary(
    admix_cov,
    af, np.array(n_samples),
    include_diag=True, abs=False
)


tile_idxs = ac.ts.create_tile_idxs(ts, tile_size=int(5e5))
sizes = [x.size for x in tile_idxs] # Number of SNPs in tiles 


straps = ac.bootstrap(
    tile_idxs, af, np.array([n_samples]*af.shape[1]).T,
    Q, ref_af, np.array([[d['n'] for d in refs]]*af.shape[1]).T, alphas,
    N_bootstrap=1e4,
    full_output = True,
    bias=True, drift_err=True,
)


out_dict = {
    'times': times,
    'Q': Q,
    'covmat': covmat,
    'admix_cov': admix_cov,
    'drift_err': drift_err,
    'G': G,
    'Ap': Ap,
    'straps_cov': straps['cov'],
    'straps_corr_cov': straps['corr_cov'],
    'straps_G_nc': straps['G_nc'],
    'straps_G': straps['G'],
    'straps_Ap': straps['Ap'],
}
with open(snakemake.output['pickle'], 'wb') as fw:
    pickle.dump(out_dict, fw)

with open(snakemake.output['text'], 'w') as fw:
    fw.writelines([
        f'Contribution of linked selection to the total variance: {G}\n',
        f'Estimated contribution of admixture to the total variance: {Ap}\n',
        "\ntile stats:\n",
        f'mean: {np.mean(sizes)}\n',
        f'std: {np.std(sizes)}\n',
        f'min: {np.min(sizes)}\n',
        f'max: {np.max(sizes)}\n',
        f'len: {len(tile_idxs)}\n',
    ])
