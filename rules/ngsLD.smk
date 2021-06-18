
rule prepare_genoFile:
  input:
    'angsd/angsd_LC_GL2_cutoff.nInd55.beagle.gz'
  output:
    'ngsLD/angsd_LC_GL2_cutoff.nInd55.geno.beagle.gz'
  log: 'log/prepare_genoFile.log'
  threads: 12
  message:
    """ Prepare beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat {input} | cut -f 4- | awk 'NR != 1' | gzip  > {output} 2> {log}
    """

rule prepare_posFile:
  input:
    'angsd/angsd_LC_GL2_cutoff.nInd55.mafs.gz'
  output:
    'ngsLD/angsd_LC_GL2_cutoff.nInd55.pos.gz'
  log: 'log/prepare_posFile.log'
  threads: 12
  message:
    """ Prepare beagle formatted genotype likelihood file generated from ANGSD (-doGlf 2) by removing the header row and the first three columns (i.e. positions, major allele, minor allele) """
  shell:
    """
    zcat {input}| cut -f 1,2 |  awk 'NR != 1' | sed 's/:/_/g'| gzip > {output}
    """

rule run_ngsLD:
  input:
   position = 'ngsLD/angsd_LC_GL2_cutoff.nInd55.pos.gz',
   geno = 'ngsLD/angsd_LC_GL2_cutoff.nInd55.geno.beagle.gz'
  output:
    'ngsLD/angsd_LC_GL2_cutoff.nInd55.ld.gz'
  log: 'log/run_ngsLD.log'
  threads: 24
  message:
    """ Estimate LD using ngsLD, which employs a maximum likelihood approach to account for genotype uncertainty in low-coverage whole genome sequencing data. """
  shell:
    """
    module load ngsld/1.1.1
    ngsld --geno {input.geno} --pos {input.position} --probs --n_ind 55 --n_sites 7239553 --max_kb_dist 0 --max_snp_dist 0 --n_threads {threads} | gzip --best > {output} 2> {log}
    """

rule run_LDpruning:
  input:
   'ngsLD/angsd_LC_GL2_cutoff.nInd55.ld.gz'
  output:
    'ngsLD/angsd_LC_GL2_cutoff.nInd55.unlinked'
  log: 'log/LD_blocks.log'
  message:
    """ Prune your data, remove SNPs with high LD """
  shell:
    """
    zcat {input} | perl scripts/prune_graph.pl --max_kb_dist 2000 --min_weight 0.5 --out {output}
    """
