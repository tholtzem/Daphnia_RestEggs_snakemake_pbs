include: "rules/common.smk"



# -----------------------------------------------

rule all:
	input:
		#expand('bbmap/rapid/{sample}.{ext}', sample=sample_names, ext=['bam']),
		##expand('getRefClones.done'),
		#expand('bbmap/rapid/minq20/{sample}.{ext}', sample=new_sample_names, ext=['minq20.bam']),
		#expand('deDup/{sample}.{ext}', sample=new_sample_names, ext=['minq20.dedup.bam', 'minq20.dedup.metrics.txt']),
		#expand('markDup/{sample}.{ext}', sample=new_sample_names, ext=['minq20.mrkdup.bam', 'minq20.mrkdup.metrics.txt']),
		#expand('deDup/{sample}.{ext}', sample=new_sample_names, ext=['minq20.overlapclipped.bam']),
		#expand('deDup/{sample}.{ext}', sample=new_sample_names, ext=['minq20.overlapclipped.bam.bai']),
		#expand('ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.fasta.{ext}', ext=['fai']),
		#expand('ref/dgal_ra_pb-target-and-other_ill-confilter.blobfilter.rmmt_sspace-lr3_lrgc3_pg3_pilon_3.simple.{ext}', ext=['dict']),
		#expand('list/overlapclippedBAM.list'),
		#expand('list/all_samples.{ext}', ext=['intervals']),
		#expand('realigned/{sample}.{ext}', sample=new_sample_names, ext=['minq20.realigned.bam']),
		#expand('depth/{sample}.{ext}', sample=new_sample_names, ext=['minq20.realigned.bam.depth.gz']),
		#expand('bedtools/{sample}.{ext}', sample=new_sample_names, ext=['minq20.realigned.genomecov.bed']),
		#expand('bedtools/plots/{sample}.{ext}', sample=new_sample_names, ext=['minq20.realigned.genomecov.pdf']),
		#expand('genome_stats.done'),
		expand('angsd/GL1_done'),
                expand('angsd/GL2_done'),
                expand('angsd/GL2_cutoffs_{IND}.done', IND=IND),
                expand('pcangsd/PCAngsd_GL1_done'),
                expand('pcangsd/PCAngsd_GL2_{IND}.done', IND=IND),
                expand('pcangsd/PCAngsd_GL2_{IND}.pdf', IND=IND),
                expand('ngsLD/angsd_LC_GL2_cutoff.nInd55.geno.beagle.gz'),
                expand('ngsLD/angsd_LC_GL2_cutoff.nInd55.pos.gz'),
                expand('ngsLD/angsd_LC_GL2_cutoff.nInd55.ld.gz')#,
                #expand('ngsLD/angsd_LC_GL2_cutoff.nInd55.unlinked')
                
		
# -----------------------------------------------


#include: "rules/s2b.smk"
include: "rules/deDup.smk"
include: "rules/indexRef.smk"
include: "rules/real_indel.smk"
include: "rules/samtools_depth.smk"
include: "rules/read_depth.smk"
include: "rules/gencov.smk"
include: "rules/angsd.smk"
include: "rules/PCAngsd.smk"
include: "rules/ngsLD.smk"
#include: "rules/angsd_cutoffs.smk"
#include: "rules/test_depth.smk"
