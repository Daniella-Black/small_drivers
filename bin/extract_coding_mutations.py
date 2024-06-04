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
from dataclasses import dataclass



@dataclass
class SequenceRange:
    """Class for the start and end of a range."""
    name: str
    start: int
    end: int
    chrom: str
    def overlaps(self, other: "SequenceRange") -> bool:
        if self.chrom != other.chrom:
            return False
        return (other.start <= self.start <= other.end) or (other.start <= self.end <= other.end) or (self.start <= other.start <= self.end) or (self.start <= other.end <= self.end)

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
regions = pd.read_csv(regions)
regions['start'] = regions['start'] -5
regions['end'] = regions['end'] +5

##chromosome notations have to be the same between regions and mutations file for overlap function so add some cleaning up in case needed (also done for mutations file)
regions['chr']  = regions['chr'].str.replace('chr','')
regions['chr']  = regions['chr'].str.replace('23','X')
regions['chr']  = regions['chr'].str.replace('23','Y')

## prepare the mutations file
mutations = pd.read_csv(mutations)
mutations['chr'] = mutations['chr'].astype('str')
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
mutations[['chr','position','REF','ALT','VAF','ASMD','CLPM','QUAL','REP']].to_csv(sample + '_coding_mutations.csv', index=False)        
