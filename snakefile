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
                #expand('angsd/angsd_GL2_cutoffs_maf0018_nInd55.done'),
                #expand('pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55.done'),
                #expand('pcangsd/PCAngsd_GL2_cutoffs_maf0018_nInd55_covmat.pdf'),
		#expand('ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_geno.beagle.gz'),
                #expand('ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub_pos.gz'),
                #expand('ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.ld.gz'),
                #expand('ngsLD/angsd_LC_GL2_cutoff_maf0018_nInd55_sub.unlinked.id'),
                #expand('LDpruned_SNPlist.done'),
                #expand('list/index_SNP.done'),
                #expand('list/getChrom_from_sites.done'),
                #expand('angsd/angsd_LC_GL2_maf0018_LDpruned.done'),
                #expand('pcangsd/PCAngsd_GL2_maf0018_LDpruned.done'),
                #expand('pcangsd/PCAngsd_GL2_maf0018_LDpruned.pdf'),
                #expand('pcangsd/PCAngsd_GL2_maf0018_LDpruned_admix_K3.done'),
                expand('angsd/{sample}.GL2.saf.idx.done', sample=BAM),
                expand('angsd/{sample}.GL2.est.ml', sample=BAM)
                


                



                
	
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
include: "rules/angsd_LDprunedSNPs.smk"
include: "rules/admix.smk"
include: "rules/realSFS.smk"
#include: "rules/angsd_cutoffs.smk"
#include: "rules/test_depth.smk"
