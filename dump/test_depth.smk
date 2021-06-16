rule angsd_depth_all:
  input:
    bamlist = 'test/data/test.list',
    ref = config["ref"]
  output:
    'test/angsd/all'
  log: 'test/log/all.log'
  threads: 12
  message:
    """ Calculate depth in ANGSD using no quality filters """
  shell:
    """
    angsd -b {input.bamlist} -doDepth 1 -maxDepth 8000 -dumpCounts 1 -out {output} -doCounts 1 -ref {input.ref} 2> {log}
    """

rule angsd_depth_strict:
  input:
    bamlist = 'test/data/test.list',
    ref = config["ref"]
  output:
    'test/angsd/strict'
  log: 'test/log/strict.log'
  threads: 12
  message:
    """ Calculate depth in ANGSD with mapping quality and nucleotide qscore above 20 """
  shell:
    """
    angsd -b {input.bamlist} -doDepth 1 -maxDepth 8000 -dumpCounts 1 -out {output} -doCounts 1 -ref {input.ref} -minMapQ 20 -minQ 20 2> {log}
    """

rule samtools_depth_all:
  input:
    'test/data/{sample}'
  output:
    'test/samtools/{sample}.depth.gz'
  log: 'log/{sample}.log'
  threads: 12
  message:
    """ Count per position depth per sample using samtools depth """
  shell:
    """
    samtools depth -aa {input} | cut -f3 | gzip > {output} 2> {log}
    """

rule samtools_depth_strict:
 input:
   'test/data/{sample}' 
 output:
   'test/samtools/{sample}.depthQ20q20.gz'
 threads: 12
 message:
    """ Count per position depth per sample using samtools depth only count reads with base quality and mapping quality greater than 20 """
 shell:
    """
    samtools depth -aa {input} -q 20 -Q 20 | cut -f3 | gzip > {output}
    """

