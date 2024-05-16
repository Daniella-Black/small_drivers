#!/usr/local/bin/python3

#####

#script to pull out all mutations within defined genomic regions of interest
#intended to pull out mutations overlapping exons
#adds/takes away 5 to each side of the exon to capture splice region variants 


######

#load packages
import pandas as pd
import re
import os
import argparse

##read in arguments
my_parser = argparse.ArgumentParser(description='get arguments')
my_parser.add_argument('-sample',
                       type=str,
                       help='sample')
my_parser.add_argument('-mutations',
                       type=str,
                       help='path to file with columns chr position REF and ALT columns for mutations (can be snv or indel)')
my_parser.add_argument('-regions',
                       type=str,
                       help='path to file with columns chr, start and end of regions of interest (e.g. exons)')

args = my_parser.parse_args()


sample = args.sample
mutations = args.mutations
regions = args.regions



###prepare the regions file: here adding 5 on either end of the exon to pick up any splice reigon variants of interest (i.e. up to +5G)
regions = pd.read_csv('/mnt/session_data/Extract_mutations_in_exons/2024_05_16_mane_and_canonical_transcript_exons_hg37_without_extra_info.csv')
regions['start'] = regions['start'] -5
regions['end'] = regions['end'] +5

##chromosome notations have to be the same between regions and mutations file for overlap function so add some cleaning up in case needed (also done for mutations file)
regions['chr']  = regions['chr'].str.replace('chr','')
regions['chr']  = regions['chr'].str.replace('23','X')
regions['chr']  = regions['chr'].str.replace('23','Y')

## prepare the mutations file
mutations = pd.read_csv('/mnt/session_data/Extract_mutations_in_exons/breast560/filtered_snv_indels/PD10010a_filtered_snv_indel_from_rdata.csv')
mutations['id'] = mutations['chr'].astype('str') + '_' +mutations['position'].astype('str') + '_'+ mutations['REF']+ '_'+ mutations['ALT']
mutations['chr']  = mutations['chr'].str.replace('chr','')
mutations['chr']  = mutations['chr'].str.replace('23','X')
mutations['chr']  = mutations['chr'].str.replace('23','Y')


##find the mutations which overlap with the regions inputted
coding = []
for row in range(len(mutations.index)):
      mutation  = SequenceRange(mutations['id'][row],  mutations['position'][row],  mutations['position'][row], mutations['chr'][row])
      regions_chr = regions[regions['chr']==mutation.chrom].reset_index(drop=True)
      for row_region in range(len(regions_chr.index)):
            region = SequenceRange('placeholder',  regions_chr['start'][row_region],  regions_chr['end'][row_region], regions_chr['chr'][row_region])
            if region.overlaps(mutation):
                coding.append(mutation.name)

##keep only the mutations overlapping with regions
mutations= mutations[mutations['id'].isin(coding)]

##output table of coding mutations
mutations[['chr','position','REF','ALT']].to_csv(sample + '_coding_mutations.csv', index=False)        
