
rule angsd_GL1:
  input:
    ref = config["ref"]
  output:
    touch('angsd/GL1_done')
  log: 'log/GL1_done.log'
  threads: 24
  message:
    """ Calculate genotype likelihoods (output: beagle) for bam files having a read depth gretaer than 10 """
  shell:
    """
    angsd -b 'list/realignedBAM_depth10.list' -ref {input.ref} -out angsd/angsd_LC_GL1 -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 20 -minQ 20 -setMinDepth 400 -setMaxDepth 2000 -doCounts 1 2> {log}
    """

rule angsd_GL2:
  input:
    ref = config["ref"]
  output:
    touch('angsd/GL2_cutoffs_{IND}.done')
  log: 'log/GL_cutoffs_{IND}.log'
  threads: 24
  message:
    """ Calculate genotype likelihoods for bam files having a read depth gretaer than 10 and try varying the cutoff -nInd and record how many sites are predicted to be variable for each scenario """
  shell:
    """
    angsd -b 'list/realignedBAM_depth10.list' -ref {input.ref} -out angsd/angsd_LC_GL2_cutoff.nInd{wildcards.IND} -GL 2 -doGlf 2 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -minMapQ 20 -minQ 20 -minInd {wildcards.IND} -setMinDepth 400 -setMaxDepth 2000 -doCounts 1 2> {log}
    """
