#!/usr/bin/env python

import pickle
import mushi
import numpy as np

# Load ksfs and true histories
ksfs = pickle.load(open('ksfs.pkl', 'rb'))
mu0 = pickle.load(open('mu0.pkl', 'rb'))

convergence_params = dict(tol=0, max_iter=${params.max_iter})
trend_kwargs = dict(max_iter=${params.trend_max_iter})

alpha_trend = ((0, ${alpha_0}), (1, ${alpha_1}))
beta_trend = ((0, ${beta_0}), (1, ${beta_1}))

ksfs.infer_eta(mu0,
               *alpha_trend,
               ridge_penalty=${alpha_ridge},
               loss='prf',
               folded=${folded},
               pts=${params.pts}, ta=${params.ta},
               trend_kwargs=trend_kwargs,
               **convergence_params, verbose=True)

ksfs.infer_mush(*beta_trend,
                ridge_penalty=${beta_ridge},
                loss='prf',
                trend_kwargs=trend_kwargs,
                **convergence_params, verbose=True)

pickle.dump([alpha_trend, beta_trend, ksfs], open('dat.pkl', 'wb'))
