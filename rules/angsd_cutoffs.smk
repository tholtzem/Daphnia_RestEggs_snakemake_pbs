rule angsd_cutoffs:
  input:
    fai = config["fai"],
    glf = 'angsd/angsd_LC_GL2.glf.gz'    
  output:
    touch('angsd/GL_cutoffs_{IND}.done')
  log: 'log/GL_cutoffs_{IND}.log'
  threads: 12
  message:
    """ Try varying the cutoff -nInd and record how many sites are predicted to be variable for each scenario.  """
  shell:
    """
    angsd -glf {input.glf} -nInd {wildcards.IND} -fai {input.fai} -out angsd/angsd_LC_GL2_cutoff.nInd{wildcards.IND} -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -skipTriallelic 1 2> {log}
    """
#echo $IND 'zcat angsd/angsd_LC_GL2_cutoff.nInd$IND.mafs.gz | tail -n+2 | wc -l' > test.txt 
