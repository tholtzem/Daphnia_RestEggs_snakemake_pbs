
rule prepare_genoFile:
  input:
    'angsd/angsd_LC_GL1.beagle.gz'
  output:
    'ngsLD/angsd_LC_GL1.fmt.beagle.gz'
  log: 'log/prepare_genoFile.log'
  threads: 12
  message:
    """ Prepare beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat {input} | cut -f 4- | gzip  > {output} 2> {log}
    """

rule prepare_posFile:
  input:
    'angsd/angsd_LC_GL1.mafs.gz'
  output:
    'ngsLD/angsd_LC_GL1.pos.gz'
  log: 'log/prepare_posFile.log'
  threads: 12
  message:
    """ Prepare beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat {input}| cut -f 1,2 | sed 's/:/_/g'| gzip > {output}
    """
