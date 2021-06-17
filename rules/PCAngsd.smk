rule PCAngsd:
  input:
    'angsd/angsd_LC_GL1.beagle.gz'
  output:
    #'angsd/PCAngsd_covmat'
    touch('angsd/PCAngsd_GL1_done')
  log: 'log/PCAngsd1_covmat.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle {input} -o pcangsd/PCAngsd_GL1_covmat -threads {threads} 2> {log}
    """

rule PCAngsd2:
  input:
    'angsd/angsd_LC_GL2_cutoff.nInd{IND}.beagle.gz'
  output:
    touch('pcangsd/PCAngsd_GL2_{IND}.done')
  log: 'log/PCAngsd_GL2_{IND}_covmat.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle {input} -o pcangsd/PCAngsd_GL2_{wildcards.IND}_covmat 2> {log}
    """

rule plot_covMat:
  input:
    touched = 'pcangsd/PCAngsd_GL2_{IND}.done',
    covMat ='pcangsd/PCAngsd_GL2_{IND}_covmat.cov'
  output:
    pdf = 'pcangsd/PCAngsd_GL2_{IND}.pdf'
  log: 'log/PCAngsd_GL2_plotCovMat{IND}.log'
  threads: 12
  message:
    """ Estimate covariance matrix from GL using PCAngsd """
  shell:
    """
    Rscript scripts/plot_covMat.R {input.covMat} {output.pdf} 2> {log}
    """

