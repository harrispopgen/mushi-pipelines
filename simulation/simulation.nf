#!/usr/bin/env nextflow

params.outdir = "output"

params.pts = 100
params.ta = 100000
params.n_haplotypes = 200
params.max_iter = 300
params.r = 0.01
params.trend_max_iter = 20

process histories {

  executor 'sge'
  memory '1 GB'
  time '1m'
  // scratch true
  conda "${CONDA_PREFIX}/envs/simulation"
  publishDir "$params.outdir", mode: 'copy'

  output:
  tuple 'eta.pkl', 'mu.pkl' into histories_ch
  tuple 'eta.pdf', 'mu.pdf' into histories_summary_ch

  """
  #!/usr/bin/env python

  import mushi
  import stdpopsim
  import numpy as np
  from matplotlib import pyplot as plt
  from scipy.special import expit
  import pickle

  # Time grid

  change_points = np.logspace(0, np.log10(${params.ta}), ${params.pts})
  t = np.concatenate((np.array([0]), change_points))


  # Sawtooth demographic history eta(t)

  species = stdpopsim.get_species("HomSap")
  model = species.get_demographic_model("Zigzag_1S14")
  ddb = model.get_demography_debugger()
  eta = mushi.eta(change_points, 2 * ddb.population_size_trajectory(steps=t)[:, 0])

  plt.figure(figsize=(3.5, 3.5))
  eta.plot(c='k')
  plt.savefig('eta.pdf')

  pickle.dump(eta, open('eta.pkl', 'wb'))


  # Mutation rate history mu(t)
  # A 96 dimensional history with a mixture of two latent signature: constant and pulse.

  flat = np.ones_like(t)
  pulse = expit(.1 * (t - 100)) - expit(.01 * (t - 2000))
  ramp = expit(-.01 * (t - 100))
  n_mutation_types = 96
  Z = np.zeros((len(t), n_mutation_types))
  np.random.seed(0)

  Z[:, 0] = 1 * flat + .5 * pulse
  Z[:, 1] = 0.5 * flat + .4 * ramp

  for col in range(2, n_mutation_types):
      scale = np.random.lognormal(-.2, .3)
      pulse_weight = 5 if col == 0 else 0
      Z[:, col] = scale * flat

  mu = mushi.mu(change_points, Z)

  plt.figure(figsize=(4, 4))
  mu.plot(range(2, n_mutation_types), alpha=0.1, lw=2, c='C0', clr=True)
  mu.plot((0,), alpha=0.75, lw=3, c='C1', clr=True)
  mu.plot((1,), alpha=0.75, lw=3, c='C2', clr=True)
  plt.savefig('mu.pdf')

  pickle.dump(mu, open('mu.pkl', 'wb'))
  """
}

process ksfs {

  executor 'sge'
  memory '10 GB'
  time '1d'
  // scratch true
  conda "${CONDA_PREFIX}/envs/simulation"
  publishDir "$params.outdir", pattern: '{*.pdf,*.txt}', mode: 'copy'

  input:
  tuple 'eta.pkl', 'mu.pkl' from histories_ch

  output:
  file 'ksfs.pkl' into ksfs_ch
  tuple 'sfs.pdf', 'ksfs.pdf', 'n_trees.txt', 'n_variants.txt', 'tmrca_cdf.pdf', 'mu0.txt' into ksfs_summary_ch
  file 'mu0.pkl' into mu0_ch

  """
  #!/usr/bin/env python

  import pickle
  import mushi
  import stdpopsim
  import msprime
  import numpy as np
  from matplotlib import pyplot as plt
  from numpy.random import binomial

  # Load histories
  eta = pickle.load(open('eta.pkl', 'rb'))
  mu = pickle.load(open('mu.pkl', 'rb'))

  # stdpopsim simulation
  # - n sampled haplotypes
  # - generate initial tree sequence without mutations, then we'll place mutations according to the MuSH

  n = ${params.n_haplotypes}
  chrom = 1
  species = stdpopsim.get_species('HomSap')
  contig = species.get_contig(f'chr{chrom}')
  contig = stdpopsim.Contig(recombination_map=contig.recombination_map, mutation_rate=0, genetic_map=contig.genetic_map)
  model = stdpopsim.PiecewiseConstantSize(eta.y[0] / 2, *zip(eta.change_points, eta.y[1:] / 2))
  samples = model.get_samples(n)
  engine = stdpopsim.get_default_engine()
  ts = engine.simulate(model, contig, samples, mutation_rate=0, seed=1)


  # Number of marginal trees
  print(ts.num_trees, file=open('n_trees.txt', 'w'))

  # Simulate k-SFS
  # - place mutations on simulated tree sequence according to mu(t)
  # - iterate over mutation types and epochs
  # - compute component of k-SFS for each iterate

  X = np.zeros((n - 1, mu.Z.shape[1]), dtype=int)
  for start_time, end_time, mutation_rate in mu.epochs():
      mutation_rate_total = mutation_rate.sum()
      print(f'epoch boundaries: ({start_time:.2f}, {end_time:.2f}), μ: {mutation_rate_total:.2f}     ',
            flush=True, end='     \\r')
      # note: the factor of 1 / chrom.length is needed to scale the mutation rate from sites to genomes
      ts_temp = msprime.mutate(ts,
                               rate=mutation_rate_total / species.genome.chromosomes[chrom - 1].length,
                               start_time=start_time,
                               end_time=end_time,
                               random_seed=1,
                               keep=True)
      for var in ts_temp.variants():
          i = var.genotypes.sum() - 1
          j = np.random.choice(mu.Z.shape[1], p=(mutation_rate / mutation_rate_total))
          X[i, j] += 1

  X_misid = binomial(X, ${params.r})
  # assume first two mutation types are misidentification partners
  X[:, 0] = X[:, 0] - X_misid[:, 0] + X_misid[::-1, 1]
  X[:, 1] = X[:, 1] - X_misid[:, 1] + X_misid[::-1, 0]
  X[:, 2:] = X[:, 2:] - X_misid[:, 2:] + X_misid[::-1, 2:]

  ksfs = mushi.kSFS(X=X)


  # Plot SFS and k-SFS

  plt.figure(figsize=(4, 3))
  ksfs.plot_total(kwargs=dict(ls='', marker='o', ms=5, c='k', alpha=0.75),
                  line_kwargs=dict(c='C0', alpha=0.75, lw=3),
                  fill_kwargs=dict(color='C0', alpha=0.1))
  plt.xscale('log')
  plt.yscale('log')
  plt.savefig('sfs.pdf')

  plt.figure(figsize=(4, 3))
  ksfs.plot(range(2, mu.Z.shape[1]), clr=True, kwargs=dict(alpha=0.1, ls='', marker='.', c='C0'))
  ksfs.plot((0,), clr=True, kwargs=dict(alpha=0.75, ls='', marker='o', c='C1'))
  ksfs.plot((1,), clr=True, kwargs=dict(alpha=0.75, ls='', marker='o', c='C2'))
  plt.savefig('ksfs.pdf')

  pickle.dump(ksfs, open('ksfs.pkl', 'wb'))

  # Number of segregating sites
  print(ksfs.X.sum(), file=open('n_variants.txt', 'w'))

  # TMRCA CDF
  plt.figure(figsize=(3.5, 3.5))
  plt.plot(eta.change_points, ksfs.tmrca_cdf(eta))
  plt.xlabel('\$t\$')
  plt.ylabel('TMRCA CDF')
  plt.ylim([0, 1])
  plt.xscale('log')
  plt.tight_layout()
  plt.savefig('tmrca_cdf.pdf')

  # Estimate constant total mutation rate using most recent time point (introducing a misspecification)
  mu0 = mu.Z[0, :].sum()
  print(mu0, file=open('mu0.txt', 'w'))
  pickle.dump(mu0, file=open('mu0.pkl', 'wb'))
  """
}

alpha_0 = [0] + (0..3.5).by(0.5).collect{ 10**it }
alpha_1 = [0] + (0..4).by(0.5).collect{ 10**it }
alpha_ridge = 1e-4

beta_0 = 1e2
beta_1 = 1e2
beta_ridge = 1e-4

folded = 'False'
freq_mask = 'False'

process eta_sweep {

  executor 'sge'
  memory '1 GB'
  time '1h'
  scratch true
  conda "${CONDA_PREFIX}/envs/simulation"
  publishDir "$params.outdir/eta_sweep/${alpha_0}_${alpha_1}", mode: 'copy'

  input:
  file 'ksfs.pkl' from ksfs_ch
  file 'mu0.pkl' from mu0_ch
  each alpha_0 from alpha_0
  each alpha_1 from alpha_1
  val alpha_ridge
  val beta_0
  val beta_1
  val beta_ridge
  val folded
  val freq_mask

  output:
  file 'dat.pkl' into eta_sweep_ch

  script:
  template 'infer.py'
}

alpha_0 = 1e2
alpha_1 = 1e3
alpha_ridge = 1e-4

beta_0 = [0] + (0..2.5).by(0.5).collect{ 10**it }
beta_1 = [0] + (0..3).by(0.5).collect{ 10**it }
beta_ridge = 1e-4

process mu_sweep {

  executor 'sge'
  memory '1 GB'
  time '1h'
  scratch true
  conda "${CONDA_PREFIX}/envs/simulation"
  publishDir "$params.outdir/mu_sweep/${beta_0}_${beta_1}", mode: 'copy'

  input:
  file 'ksfs.pkl' from ksfs_ch
  file 'mu0.pkl' from mu0_ch
  val alpha_0
  val alpha_1
  val alpha_ridge
  each beta_0 from beta_0
  each beta_1 from beta_1
  val beta_ridge
  val folded
  val freq_mask

  output:
  file 'dat.pkl' into mu_sweep_ch

  script:
  template 'infer.py'
}


beta_0 = 1e2
beta_1 = 1e2
beta_ridge = 1e-4

folded = ['False', 'True']
freq_mask = 'True'

process mush {

  executor 'sge'
  memory '1 GB'
  time '1h'
  scratch true
  conda "${CONDA_PREFIX}/envs/simulation"
  publishDir "$params.outdir/mush/folded_${folded}", mode: 'copy'

  input:
  file 'ksfs.pkl' from ksfs_ch
  each folded from folded
  file 'mu0.pkl' from mu0_ch
  val alpha_0
  val alpha_1
  val alpha_ridge
  val beta_0
  val beta_1
  val beta_ridge
  val freq_mask

  output:
  file 'dat.pkl' into mush_ch

  script:
  template 'infer.py'
}
