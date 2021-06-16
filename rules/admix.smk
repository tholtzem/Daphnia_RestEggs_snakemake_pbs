rule PCAngsd_admix:
  input:
    #'angsd/angsd.beagle.gz'
  output:
    'angsd/PCAngsd'
  log: 'log/PCAngsd.log'
  threads: 12
  message:
    """ Infer admixture proportions using PCAngsd """
  shell:
    """
    module load pcangsd/1.01
    pcangsd -beagle angsd/angsd.beagle.gz -admix -o {output} -threads {threads}
    """

