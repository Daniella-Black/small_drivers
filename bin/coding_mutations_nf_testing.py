#!/usr/local/bin/python3
import pandas as pd
import re
import os
import argparse

my_parser = argparse.ArgumentParser(description='get arguments')
my_parser.add_argument('-sample',
                       type=str,
                       help='sample')
my_parser.add_argument('-annotation_vcf_path',
                       type=str,
                       help='annotation_vcf_path')
my_parser.add_argument('-mane',
                       type=str,
                       help='path to the mane transcript table')
my_parser.add_argument('-cmc',
                       type=str,
                       help='path to the cmc table')
my_parser.add_argument('-hgnc',
                       type=str,
                       help='path to the hgnc table')
my_parser.add_argument('-non_mane_transcripts',
                       type=str,
                       help='path to table containing non-mane transcript info')
my_parser.add_argument('-cgc',
                       type=str,
                       help='path to cgc table')
args = my_parser.parse_args()

###pull out coding mutations
sample = args.sample
annotation_vcf_path = args.annotation_vcf_path
mane_path = args.mane
cmc = args.cmc
hgnc = args.hgnc
non_mane_transcripts = args.non_mane_transcripts
cgc = args.cgc

def flatten(A):
    rt = []
    for i in A:
        if isinstance(i,list): rt.extend(flatten(i))
        else: rt.append(i)
    return rt

##terms of interested from annoated vcf column for drivers 
relevant_terms = ['splice_acceptor_variant',
'splice_donor_variant',
'stop_gained',
'frameshift_variant',
'stop_lost',
'start_lost',
'inframe_insertion',
'inframe_deletion',
'missense_variant',
#'splice_donor_5th_base_variant'
'splice_region_variant']

##read in the cgc table and get the ensg numbers into correct format
cgc = pd.read_csv(cgc, sep='\t')
#get the cleaned up ensg number                       
cgc[['start', 'containing_ensg']] = cgc.Synonyms.str.split('ENSG', expand=True)
cgc[['ensg_no', 'rest']] = cgc.containing_ensg.str.split(".", 1, expand=True)
cgc['ENSG'] = 'ENSG' + cgc['ensg_no']
##tidy up
cgc = cgc.drop(['start', 'containing_ensg','ensg_no', 'rest'], axis=1)

#get the mane file into the correct format for nextflow
mane = pd.read_csv(mane_path, sep='\t')
mane[['transcript_ID', 'to_del2']] = mane['name'].str.split('.', 1,expand=True)
mane[['gene_ID', 'to_del2']] = mane['geneName'].str.split('.', 1, expand=True)
mane = mane.rename(columns={'#chrom': 'chr', 'chromStart': 'start', 'chromEnd': 'end', 'geneName2': 'gene_name' })
mane['chr'] = mane['chr'].str.replace('chr', '')
mane_full = mane                       
mane = mane[['chr', 'start', 'end', 'transcript_ID','gene_ID', 'gene_name']]
mane.to_csv('mane_transcripts.csv',index=False)

non_mane_transcripts = pd.read_csv(non_mane_transcripts,sep='\t')
non_mane_transcripts.to_csv('non_mane_transcripts.csv',index=False)
print(str(len(mane.index)))
print(str(len(non_mane_transcripts.index)))
mane = mane[~mane['gene_ID'].isin(list(non_mane_transcripts['gene_ID']))]
print(str(len(mane.index)))
mane = pd.concat([mane, non_mane_transcripts])
print(str(len(mane.index)))
mane = mane.reset_index(drop = True)

mane.to_csv('merge.csv',index=False)
