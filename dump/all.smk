rule mergeBAM:
  input:
    aln = 'bbmap/rapid/{ID}.bb.RAPID.bam'
  output:
    mrgd = 'bbmap/rapid/{ID}.bb.RAPID_merged.bam'
  threads: 12
  message: """--- Merging bam files with samtools merge ---"""
  shell:
    """
    samtools merge {output.mrgd} {input.aln}
    """

rule mark_duplicates:
  input:
    mrgd = mrgd = 'bbmap/rapid/{ID}.bb.RAPID_merged.bam'
  output:
    mrkdDup = 'bbmap/rapid/{ID}.bb.RAPID.mrkddups.bam'
  threads: 12
  message: """--- Marking duplicates of merged bam files with Picard ---"""
  shell:
      """
      MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.mrgd} OUTPUT={output.mrkdDup} METRICS_FILE={output.mrkdDup}.metrics
      """

rule remove_duplicates:
  input:
    mrgd = mrgd = 'bbmap/rapid/{ID}.bb.RAPID_merged.bam'
  output:
    rmDup = 'bbmap/rapid/{ID}.bb.RAPID.dedup.bam'
  threads: 12
  message: """--- Removing duplicates of merged bam files with Picard ---"""
  shell:
      """
      MarkDuplicates TMP_DIR=./ MAX_RECORDS_IN_RAM=15000000 REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 CREATE_INDEX=true INPUT={input.mrgd} OUTPUT={output.rmDup} METRICS_FILE={output.rmDup}.metrics
      """


rule bamIndex:
  input:
    rmDup = 'bbmap/rapid/{ID}.bb.RAPID.dedup.bam'
  output:
    idx = 'bbmap/rapid/{ID}.bb.RAPID.dedup.bam.bai'
  threads: 2
  message: """--- Indexing with samtools ---"""
  shell:
    """
    samtools index {input.aln} {output.idx}
    """
