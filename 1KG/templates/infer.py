#!/usr/bin/env python

import pickle
import mushi
import pandas as pd
import numpy as np

ksfs = mushi.kSFS(file='ksfs.tsv')

# pre-specified eta
if ${eta}:
    eta = pickle.load(open('eta.pkl', 'rb'))
else:
    eta = None

# reference eta for ancestral fusion
if ${ref_pop}:
    ref_pop = pickle.load(open('dat.ref.pkl', 'rb'))[2]
    eta_ref = ref_pop.eta
    mu_ref = ref_pop.mu
else:
    eta_ref = None
    mu_ref = None

# sorts the columns of the ksfs (and other relevant attributes)
sorted_triplets = [f'{a5}{a}{a3}>{a5}{d}{a3}' for a in 'AC'
                   for d in 'ACGT' if d != a
                   for a5 in 'ACGT' for a3 in 'ACGT']
foo, bar = ksfs.mutation_types.reindex(sorted_triplets)
ksfs.mutation_types = foo
ksfs.X = ksfs.X[:, bar]
ksfs.AM_mut = ksfs.AM_mut[bar, :][:, bar]

# TCC>TTC and GAA>GGA only
if ${tcc}:
  idx = ksfs.mutation_types.get_loc('TCC>TTC')
  idx_rev = ksfs.mutation_types.get_loc('GAA>GGA')

  X = np.zeros((ksfs.X.shape[0], 3))
  X[:, 0] = ksfs.X[:, idx]
  X[:, 1] = ksfs.X[:, idx_rev]
  X[:, 2] = ksfs.X.sum(1) - X[:, 0] - X[:, 1]

  mutation_types = ['TCC>TTC', 'GAA>GGA', None]

  ksfs = mushi.kSFS(X=X, mutation_types=mutation_types)

if ${boot}:
    from scipy.stats import multinomial
    p = np.ravel(np.array(ksfs.X, dtype=float))
    p /= p.sum()
    n = ksfs.X.sum()
    ksfs.X = multinomial.rvs(n, p).reshape(ksfs.X.shape)

masked_genome_size = pd.read_csv('masked_size.tsv', sep='\t', header=None,
                                 index_col=0, names=('count',))

u = ${params.mut_rate}
mu0 = u * masked_genome_size['count'].sum()

convergence_params = dict(tol=0, max_iter=${params.max_iter})
trend_kwargs = dict(max_iter=${params.trend_max_iter})

if ${k_eta2} is None:
    alpha_trend = ((${k_eta1}, ${lambda_eta1}),)
else:
    alpha_trend = ((${k_eta1}, ${lambda_eta1}), (${k_eta2}, ${lambda_eta2}))

dat = []

if eta is None:
    ksfs.infer_eta(mu0,
                   *alpha_trend,
                   folded=${folded},
                   ridge_penalty=${alpha_ridge},
                   eta_ref=eta_ref,
                   loss='prf',
                   pts=${params.pts}, ta=${params.ta},
                   trend_kwargs=trend_kwargs,
                   **convergence_params, verbose=True)
    dat.append(alpha_trend)
else:
    ksfs.set_eta(eta)
    ksfs.infer_misid(mu0, verbose=True)

if ${folded}:
    beta_trend = None
else:
    if ${k_mu2} is None:
        beta_trend = ((${k_mu1}, ${lambda_mu1}),)
    else:
        beta_trend = ((${k_mu1}, ${lambda_mu1}), (${k_mu2}, ${lambda_mu2}))

    ksfs.infer_mush(*beta_trend,
                    ridge_penalty=${beta_ridge},
                    rank_penalty=${beta_rank},
                    mu_ref=mu_ref,
                    loss='prf',
                    trend_kwargs=trend_kwargs,
                    **convergence_params, verbose=True)

dat += [beta_trend, ksfs, '${population}']

pickle.dump(dat, open('dat.pkl', 'wb'))
