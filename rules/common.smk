import os
import pandas as pd

###### Config file and sample sheets #####


configfile: 
  "config/config.yaml"

# load sample info from sam files 
samples_information = pd.read_csv("list/samRAPID.txt", sep='\t', index_col=False)
sample_names = list(samples_information['sample'])
#sample_locations = list(samples_information['location'])
sample_ID = list(samples_information['id'])

# load sample info from bam files
new_samples_information = pd.read_csv("list/newRAPID.txt", sep='\t', index_col=False)
new_sample_names = list(new_samples_information['sample'])
new_sample_ID = list(new_samples_information['id'])


# load number of individuals for angsd_cutoffs
nInd = pd.read_csv("list/cutoff_nInd.txt", sep='\t')
IND = list(nInd['nInd'])

prefixBAMdepth10 = pd.read_csv("list/prefix_realignedBAM_depth10.list", sep='\t')
BAM = list(prefixBAMdepth10['prefix'])
