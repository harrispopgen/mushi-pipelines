#!/usr/bin/env nextflow

// params.vcf_dir = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/"
params.vcf_dir = "/net/harris/vol1/data/30x1KG/"
params.mask = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed"
params.ancestor = "ftp://ftp.ensembl.org/pub/release-100/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh38.tar.gz"
params.samples = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
// params.outdir = "output_phase3"
params.outdir = "output_30x"
params.k = 3

params.mut_rate = 1.25e-8

params.pts = 200
params.ta = 200000
params.max_iter = 500
params.trend_max_iter = 50

chromosomes = 1..22

// vcf_ch = Channel
//   .of (chromosomes)
//   .map { [it,
//           file(params.vcf_dir + "ALL.chr${it}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"),
//           file(params.vcf_dir + "ALL.chr${it}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi")] }

vcf_ch = Channel
  .of (chromosomes)
  .map { [it,
          file(params.vcf_dir + "CCDG_13607_B01_GRM_WGS_2019-02-19_chr${it}.recalibrated_variants.vcf.gz"),
          file(params.vcf_dir + "CCDG_13607_B01_GRM_WGS_2019-02-19_chr${it}.recalibrated_variants.vcf.gz.tbi")] }

process mask {

  executor 'sge'
  memory '10 MB'
  time '10m'
  scratch true

  input:
  path 'mask.allchr.bed' from params.mask
  each chromosome from chromosomes

  output:
  tuple chromosome, 'mask.bed' into mask_ch

  """
  grep -P "^chr${chromosome}\\t" mask.allchr.bed | cut -f1-3 > mask.bed
  """
}

process ancestor {

  executor 'sge'
  memory '10 MB'
  time '10m'
  scratch true

  input:
  path 'homo_sapiens_ancestor_GRCh38.tar.gz' from params.ancestor
  each chromosome from chromosomes

  output:
  tuple chromosome, 'ancestor.fa' into ancestor_ch

  """
  tar -zxvf homo_sapiens_ancestor_GRCh38.tar.gz homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${chromosome}.fa
  echo ">chr${chromosome}" > ancestor.fa
  tail -n +2 homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${chromosome}.fa >> ancestor.fa
  """
}

mask_ch.into{ mask_ch_1; mask_ch_2 }
ancestor_ch.into{ ancestor_ch_1; ancestor_ch_2 }

process masked_size {

  executor 'sge'
  memory '100 MB'
  time '10h'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"

  input:
  tuple chrom, 'mask.bed', 'ancestor.fa' from mask_ch_1.join(ancestor_ch_1)
  val k from params.k

  output:
  file 'masked_size.tsv' into masked_size_ch

  """
  mutyper targets ancestor.fa --strict --k ${k} --bed mask.bed > masked_size.tsv
  """
}

process masked_size_total {

  executor 'sge'
  memory '100 MB'
  time '10m'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir params.outdir, mode: 'copy'

  input:
  file 'masked_size' from masked_size_ch.collect()

  output:
  file 'masked_size.tsv' into masked_size_total_ch

  """
  #! /usr/bin/env python

  import glob
  import pandas as pd

  sum(pd.read_csv(file, sep='\t', index_col=0, header=None, squeeze=True)
      for file in glob.glob('masked_size*')).to_csv('masked_size.tsv', sep='\t', header=False)
  """
}

// mutation types for each chromosome vcf
process mutation_types {

  executor 'sge'
  memory '500 MB'
  time '2d'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"

  input:
  tuple chrom, 'mask.bed', 'ancestor.fa', 'snps.vcf.gz', 'snps.vcf.gz.tbi' from mask_ch_2.join(ancestor_ch_2).join(vcf_ch)
  path 'integrated_call_samples.tsv' from params.samples
  val k from params.k

  output:
  tuple chrom, 'mutation_types.vcf.gz' into mutation_types_ch

  """
  # NOTE: sample NA18498 is missing in hg38 release of low coverage data, see https://www.biorxiv.org/content/10.1101/600254v1
  # tail -n +2 integrated_call_samples.tsv | cut -f1 | grep -v NA18498 > all_samples.txt
  tail -n +2 integrated_call_samples.tsv | cut -f1 > all_samples.txt
  bcftools view -S all_samples.txt -c 1:minor -R mask.bed -m2 -M2 -v snps -f PASS -Ou snps.vcf.gz | bcftools view -g ^miss -Ou | mutyper variants ancestor.fa - --strict --k ${k} | bcftools convert -Oz -o mutation_types.vcf.gz
  """
}

populations_ch = Channel
    .fromPath(params.samples)
    .splitCsv(skip: 1, sep: '\t')
    .map{ row -> row[2] + '_' + row[1] }
    .unique()
    .into { populations_ch_1; populations_ch_2 }


process ksfs {

  executor 'sge'
  memory '500 MB'
  time '1d'
  // scratch true
  conda "${CONDA_PREFIX}/envs/1KG"

  input:
  tuple chrom, 'mutation_types.vcf.gz' from mutation_types_ch
  each pop from populations_ch_1
  path 'integrated_call_samples.tsv' from params.samples
  val k from params.k

  output:
  tuple pop, 'ksfs.tsv' into ksfs_ch

  """
  # NOTE: sample NA18498 is missing in hg38 release of low coverage data, see https://www.biorxiv.org/content/10.1101/600254v1
  # grep ${pop.split('_')[1]} integrated_call_samples.tsv | cut -f1 | grep -v NA18498 > samples.txt
  grep ${pop.split('_')[1]} integrated_call_samples.tsv | cut -f1 > samples.txt
  bcftools view -S samples.txt -c 1:minor -G mutation_types.vcf.gz | mutyper ksfs - > ksfs.tsv
  """
}

// ksfs for each population
process ksfs_total {

  executor 'sge'
  memory '100 MB'
  time '10m'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir "$params.outdir/ksfs/${pop}", mode: 'copy'

  input:
  tuple pop, 'ksfs' from ksfs_ch.groupTuple(size: chromosomes.size())

  output:
  tuple pop, 'ksfs.tsv' into ksfs_total_ch

  """
  #! /usr/bin/env python

  import glob
  import pandas as pd

  sum(pd.read_csv(file, sep='\t', index_col=0)
      for file in glob.glob('ksfs*')).to_csv('ksfs.tsv', sep='\t')
  """
}

ksfs_total_ch.into { ksfs_total_ch_1; ksfs_total_ch_2; ksfs_total_ch_3; ksfs_total_ch_4; ksfs_total_ch_5; ksfs_total_ch_6; ksfs_total_ch_7; ksfs_total_ch_8; ksfs_total_ch_9; ksfs_total_ch_10 }

// trend orders and penalty strengths
k_eta1 = [0, 1, 2, 3]
lambda_eta1 = [0] + (0..4).by(0.25).collect { 10**it }
k_eta2 = 'None'
lambda_eta2 = 'None'

alpha_ridge = 1e-4

k_mu1 = 0
lambda_mu1 = 5e1
k_mu2 = 3
lambda_mu2 = 1e-4

beta_ridge = 1e-4
beta_rank = 0

ref_pop = 'False'
folded = 'False'
eta = 'False'

boot = 'False'
tcc = 'True'

process eta_sweep {

  executor 'sge'
  memory '500 MB'
  time '1d'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir "$params.outdir/eta_sweep/${k_eta1}_${lambda_eta1}/${population}", mode: 'copy'

  input:
  tuple population, 'ksfs.tsv' from ksfs_total_ch_1.filter { it[0] == 'EUR_CEU' }
  file 'masked_size.tsv' from masked_size_total_ch
  each k_eta1 from k_eta1
  each lambda_eta1 from lambda_eta1
  val k_eta2
  val lambda_eta2
  val alpha_ridge
  val k_mu1
  val lambda_mu1
  val k_mu2
  val lambda_mu2
  val beta_ridge
  val beta_rank
  val ref_pop
  val folded
  val eta
  val boot
  val tcc

  output:
  file 'dat.pkl' into eta_sweep_ch

  script:
  template 'infer.py'
}

k_eta1 = 0
lambda_eta1 = 4e2
k_eta2 = 3
lambda_eta2 = 1e-1

k_mu1 = [0, 1, 2, 3]
lambda_mu1 = [0] + (0..4).by(0.25).collect { 10**it }

process mu_sweep {

  executor 'sge'
  memory '500 MB'
  time '1d'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir "$params.outdir/mu_sweep/${k_mu1}_${lambda_mu1}/${population}", mode: 'copy'

  input:
  tuple population, 'ksfs.tsv' from ksfs_total_ch_2.filter { it[0] == 'EUR_CEU' }
  file 'masked_size.tsv' from masked_size_total_ch
  val k_eta1
  val lambda_eta1
  val k_eta2
  val lambda_eta2
  val alpha_ridge
  each k_mu1 from k_mu1
  each lambda_mu1 from lambda_mu1
  val k_mu2
  val lambda_mu2
  val beta_ridge
  val beta_rank
  val ref_pop
  val folded
  val eta
  val boot
  val tcc

  output:
  file 'dat.pkl' into mu_sweep_ch

  script:
  template 'infer.py'
}

k_mu1 = 0
lambda_mu1 = 5e1

boot = 'True'

process bootstrap {

  executor 'sge'
  memory '500 MB'
  time '1d'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir "$params.outdir/bootstrap/${bootstrap}", mode: 'copy'

  input:
  tuple population, 'ksfs.tsv' from ksfs_total_ch_3.filter { it[0] == 'EUR_CEU' }
  file 'masked_size.tsv' from masked_size_total_ch
  each bootstrap from 1..20
  val k_eta1
  val lambda_eta1
  val k_eta2
  val lambda_eta2
  val alpha_ridge
  val k_mu1
  val lambda_mu1
  val k_mu2
  val lambda_mu2
  val beta_ridge
  val beta_rank
  val ref_pop
  val folded
  val eta
  val boot
  val tcc

  output:
  file 'dat.pkl' into bootstrap_ch

  script:
  template 'infer.py'
}

boot = 'False'

process europulse {

  executor 'sge'
  memory '500 MB'
  time '1d'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir "$params.outdir/europulse_mushi/${population}", mode: 'copy'

  input:
  tuple population, 'ksfs.tsv' from ksfs_total_ch_4.filter { it[0].split('_')[0] == 'EUR' }
  file 'masked_size.tsv' from masked_size_total_ch
  val k_eta1
  val lambda_eta1
  val k_eta2
  val lambda_eta2
  val alpha_ridge
  val k_mu1
  val lambda_mu1
  val k_mu2
  val lambda_mu2
  val beta_ridge
  val beta_rank
  val ref_pop
  val folded
  val eta
  val boot
  val tcc

  output:
  file 'dat.pkl' into europulse_ch

  script:
  template 'infer.py'
}

process eta_Tennessen {
  executor 'sge'
  memory '100 MB'
  time '1h'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"

  output:
  file 'eta.pkl' into eta_Tennessen_ch

  """
  #! /usr/bin/env python

  import numpy as np
  import stdpopsim
  import mushi
  import pickle

  species = stdpopsim.get_species("HomSap")
  model = species.get_demographic_model("OutOfAfrica_2T12")
  ddb = model.get_demography_debugger()
  change_points = np.logspace(np.log10(1), np.log10(200000), 200)
  steps = np.concatenate((np.array([0]), change_points))
  eta = mushi.eta(change_points,
                  1 / ddb.coalescence_rate_trajectory(steps=steps,
                                                      num_samples=[0, 2],
                                                      double_step_validation=False)[0])
  pickle.dump(eta, open('eta.pkl', 'wb'))
  """
}

relate_files_ch = populations_ch_2.filter { it.split('_')[0] == 'EUR' }
  .map { [it,
          file("Relate_histories/relate_${it.split('_')[1]}.coal")] }

process eta_Relate {

  executor 'sge'
  memory '100 MB'
  time '1h'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"

  input:
  tuple pop, 'file.coal' from relate_files_ch

  output:
  tuple pop, 'eta.pkl' into eta_Relate_ch

  """
  #! /usr/bin/env python

  import numpy as np
  import mushi
  import scipy
  import pickle

  pop = '${pop}'.split('_')[1]

  with open('file.coal') as f:
      f.readline()
      t = np.fromstring(f.readline(), sep=' ')
      with np.errstate(divide='ignore'):
          y = 1 / np.fromstring(f.readline(), sep=' ')[2:]
  change_points = np.logspace(np.log10(1), np.log10(200000), 200)
  t2 = np.concatenate((np.array([0]), change_points, np.array([np.inf])))
  eta = mushi.eta(change_points, scipy.interpolate.interp1d(t, y, kind='nearest')(t2[:-1]))

  pickle.dump(eta, open('eta.pkl', 'wb'))
  """
}

eta = 'True'

// more penalty for Tennessen and Relate to get closer to pulse shape
k_mu1 = 0
lambda_mu1 = 1e2

process europulse_Tennessen {

  executor 'sge'
  memory '500 MB'
  time '1d'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir "$params.outdir/europulse_Tennessen/${population}", mode: 'copy'

  input:
  tuple population, 'ksfs.tsv' from ksfs_total_ch_5.filter { it[0].split('_')[0] == 'EUR' }
  file 'masked_size.tsv' from masked_size_total_ch
  file 'eta.pkl' from eta_Tennessen_ch
  val k_eta1
  val lambda_eta1
  val k_eta2
  val lambda_eta2
  val alpha_ridge
  val k_mu1
  val lambda_mu1
  val k_mu2
  val lambda_mu2
  val beta_ridge
  val beta_rank
  val ref_pop
  val folded
  val eta
  val boot
  val tcc

  output:
  file 'dat.pkl' into europulse_Tennessen_ch

  script:
  template 'infer.py'
}

ksfs_total_ch_EUR = ksfs_total_ch_6.filter { it[0].split('_')[0] == 'EUR' }

process europulse_Relate {

  executor 'sge'
  memory '500 MB'
  time '1d'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir "$params.outdir/europulse_Relate/${population}", mode: 'copy'

  input:
  tuple population, 'ksfs.tsv', 'eta.pkl' from ksfs_total_ch_EUR.join(eta_Relate_ch)
  file 'masked_size.tsv' from masked_size_total_ch
  val k_eta1
  val lambda_eta1
  val k_eta2
  val lambda_eta2
  val alpha_ridge
  val k_mu1
  val lambda_mu1
  val k_mu2
  val lambda_mu2
  val beta_ridge
  val beta_rank
  val ref_pop
  val folded
  val eta
  val boot
  val tcc

  output:
  file 'dat.pkl' into europulse_Relate_ch

  script:
  template 'infer.py'
}

eta = 'False'
tcc = 'False'

// same as above, but all populations, ancestral fusion to YRI, rank penalty, and softer mutation spectrum history
k_mu1 = 0
lambda_mu1 = 2e2
k_mu2 = 3
lambda_mu2 = 1e0

beta_rank = 1e2

process mush_ref {

  executor 'sge'
  memory '500 MB'
  time '1d'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir "$params.outdir/mush/${population}", mode: 'copy'

  input:
  tuple population, 'ksfs.tsv' from ksfs_total_ch_7.first { it[0] == 'AFR_YRI' }
  file 'masked_size.tsv' from masked_size_total_ch
  val k_eta1
  val lambda_eta1
  val k_eta2
  val lambda_eta2
  val alpha_ridge
  val k_mu1
  val lambda_mu1
  val k_mu2
  val lambda_mu2
  val beta_ridge
  val beta_rank
  val ref_pop
  val folded
  val eta
  val boot
  val tcc

  output:
  file 'dat.pkl' into mush_ref_ch

  script:
  template 'infer.py'
}

alpha_ridge = 1e5
beta_ridge = 1e2
ref_pop = 'True'
process mush {

  executor 'sge'
  memory '500 MB'
  time '1d'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir "$params.outdir/mush/${population}", mode: 'copy'

  input:
  tuple population, 'ksfs.tsv' from ksfs_total_ch_8.filter { it[0] != 'AFR_YRI' }
  file 'masked_size.tsv' from masked_size_total_ch
  val k_eta1
  val lambda_eta1
  val k_eta2
  val lambda_eta2
  val alpha_ridge
  val k_mu1
  val lambda_mu1
  val k_mu2
  val lambda_mu2
  val beta_ridge
  val beta_rank
  val ref_pop
  file 'dat.ref.pkl' from mush_ref_ch
  val folded
  val eta
  val boot
  val tcc

  output:
  file 'dat.pkl' into mush_ch

  script:
  template 'infer.py'
}

// same as previous two, but folded
alpha_ridge = 1e-4
beta_ridge = 1e-4
ref_pop = 'False'
folded = 'True'
process mush_ref_folded {

  executor 'sge'
  memory '500 MB'
  time '1d'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir "$params.outdir/mush_folded/${population}", mode: 'copy'

  input:
  tuple population, 'ksfs.tsv' from ksfs_total_ch_9.first { it[0] == 'AFR_YRI' }
  file 'masked_size.tsv' from masked_size_total_ch
  val k_eta1
  val lambda_eta1
  val k_eta2
  val lambda_eta2
  val alpha_ridge
  val k_mu1
  val lambda_mu1
  val k_mu2
  val lambda_mu2
  val beta_ridge
  val beta_rank
  val ref_pop
  val folded
  val eta
  val boot
  val tcc

  output:
  file 'dat.pkl' into mush_ref_folded_ch

  script:
  template 'infer.py'
}

alpha_ridge = 1e4
beta_ridge = 1e2
ref_pop = 'True'
process mush_folded {

  executor 'sge'
  memory '500 MB'
  time '1d'
  scratch true
  conda "${CONDA_PREFIX}/envs/1KG"
  publishDir "$params.outdir/mush_folded/${population}", mode: 'copy'

  input:
  tuple population, 'ksfs.tsv' from ksfs_total_ch_10.filter { it[0] != 'AFR_YRI' }
  file 'masked_size.tsv' from masked_size_total_ch
  val k_eta1
  val lambda_eta1
  val k_eta2
  val lambda_eta2
  val alpha_ridge
  val k_mu1
  val lambda_mu1
  val k_mu2
  val lambda_mu2
  val beta_ridge
  val beta_rank
  val ref_pop
  file 'dat.ref.pkl' from mush_ref_folded_ch
  val folded
  val eta
  val boot
  val tcc

  output:
  file 'dat.pkl' into mush_folded_ch

  script:
  template 'infer.py'
}
