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

args = my_parser.parse_args()

###pull out coding mutations
sample = args.sample
annotation_vcf_path = args.annotation_vcf_path
mane_path = args.mane
cmc = args.cmc
hgnc = args.hgnc

def flatten(A):
    rt = []
    for i in A:
        if isinstance(i,list): rt.extend(flatten(i))
        else: rt.append(i)
    return rt

##terms of interested from annoated vcf INFO column for drivers 
relevant_terms = ['splice_acceptor_variant',
'splice_donor_variant',
'stop_gained',
'frameshift_variant',
'stop_lost',
'start_lost',
'inframe_insertion',
'inframe_deletion',
'missense_variant',
'splice_donor_5th_base_variant']


def nth_repl(s, sub, repl, n):
    find = s.find(sub)
    # If find is not -1 we have found at least one match for the substring
    i = find != -1
    # loop util we find the nth or we find no match
    while find != -1 and i != n:
        # find + 1 means we start searching from after the last match
        find = s.find(sub, find + 1)
        i += 1
    # If i is equal to n we found nth match so replace
    if i == n:
        return s[:find] + repl + s[find+len(sub):]
    return s

  

#samp = pd.read_csv('/home/jovyan/session_data/mounted-data/LP3000429-DNA_G03.vcf', comment='#', sep='\t', header=None)
samp = pd.read_csv(annotation_vcf_path, comment='#', sep='\t', header=None)
samp = samp.rename(columns={0: 'chr', 1:'pos', 2:'ID', 3:'REF', 4:'ALT', 5:'QUAL', 6:'FILTER', 7:'INFO', 8:'FORMAT',9:'NORMAL', 10:'TUMOR'})
#CHROM POSIDREF ALT QUAL FILTER INFO FORMAT
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
samp = samp[samp['chr'].isin(chroms)]
samp = samp[samp['FILTER']=='PASS']
effects = list()
for variant in range(len(samp)):
  effects.append(flatten(re.findall(r'\|\|(.*?)\|\|\|\|', samp['INFO'][variant])))
varianteffectsdf = pd.DataFrame()
varianteffectsdf[sample + '_effect'] = flatten(effects)
varianteffectsdf.to_csv(sample + '_all_variant_effects.csv')                                                                       
