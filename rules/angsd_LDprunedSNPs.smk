rule angsd_index_sites:
  input:
    touched = 'LDpruned_SNPlist.done',
    snps = "list/LDpruned_snps.list"
  output:
    touch('angsd/index_SNP.done')
  log: 'log/index_SNP.done'
  threads: 12
  message:
    """ Index SNP list using angsd """
  shell:
    """
    angsd sites index {input.snps} 2> {log}
    """

rule angsd_GL2_LDpruned:
  input:
    ref = config["ref"],
    touched = 'angsd/index_SNP.done',
    snps = "list/LDpruned_snps.list"
  output:
    touch('angsd/angsd_LC_GL2_LDpruned.done')
  log: 'log/angsd_LC_GL2_LDpruned.log'
  threads: 12
  message:
    """ Calculate genotype likelihoods on LD pruned sites only using angsd """
  shell:
    """
    angsd -b 'list/realignedBAM_depth10.list' -ref {input.ref} -out angsd/angsd_LC_GL2_LDpruned -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doCounts 1 -sites {input.snps} 2> {log}
    """

rule PCAngsd_LDpruned:
  input:
    touched = 'angsd/angsd_LC_GL2_LDpruned.done',
    #GL2 = 'angsd/angsd_LC_GL2_LDpruned.beagle.gz'
  output:
    touch('pcangsd/PCAngsd_GL2_LDpruned.done')
  log: 'log/PCAngsd_GL2_LDpruned_covmat.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle angsd/angsd_LC_GL2_LDpruned.beagle.gz -o pcangsd/PCAngsd_GL2_LDpruned_covmat 2> {log}
    """

rule plot_covMat_LDpruned:
  input:
    touched = 'pcangsd/PCAngsd_GL2_LDpruned.done',
    #covMat = 'pcangsd/PCAngsd_GL2_LDpruned_covmat.cov'
  output:
    pdf = 'pcangsd/PCAngsd_GL2_LDpruned.pdf'
  log: 'log/PCAngsd_GL2_LDpruned_plotCovMat.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    Rscript scripts/plot_covMat.R pcangsd/PCAngsd_GL2_LDpruned_covmat.cov {output.pdf} 2> {log}
    """
