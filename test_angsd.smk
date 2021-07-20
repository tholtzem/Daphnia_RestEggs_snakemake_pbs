rule test_angsd935:
  input:
    ref = config["ref"],
    #touched = 'list/index_snp_dgal1.done',
    #touched2 = 'list/getChrom_dgal1_from_sites.done'
    touched = 'list/index_SNP.done',
    touched2 = 'list/getChrom_from_sites.done'
  output:
    touch('test_angsd935.done')
    #touch('angsd/angsd_LC_GL2_maf0018_LDpruned_dgal1.done')
  log:
    'test_angsd935.log'
    #'log/angsd_LC_GL2_maf0018_LDpruned_dgal1.log'
  threads: 24
  message:
    """ Calculate genotype likelihoods on predefined LD pruned sites only using angsd (-doMajorMinor 3 -sites)  """
  shell:
    """
    module load angsd/0.935
    angsd -b 'list/realignedBAM_depth10.list' -anc {input.ref} -out test_angsd935 -GL 2 -doGlf 2 -doMajorMinor 3 -doMaf 1 -doCounts 1 -sites 'list/LDpruned_snps.list' -rf list/LDpruned_Chrom.list 2> {log}
    """


