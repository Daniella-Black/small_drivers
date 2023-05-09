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
args = my_parser.parse_args()

###pull out coding mutations
sample = args.sample
annotation_vcf_path = args.annotation_vcf_path
mane_path = args.mane
cmc = args.cmc
hgnc = args.hgnc
non_mane_transcripts = args.non_mane_transcripts

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

  
cmc_colnames =['GENE_NAME', 'ACCESSION_NUMBER', 'ONC_TSG', 'CGC_TIER', 'MUTATION_URL',
       'LEGACY_MUTATION_ID', 'Mutation CDS', 'Mutation AA', 'AA_MUT_START',
       'AA_MUT_STOP', 'SHARED_AA', 'GENOMIC_WT_ALLELE_SEQ',
       'GENOMIC_MUT_ALLELE_SEQ', 'AA_WT_ALLELE_SEQ', 'AA_MUT_ALLELE_SEQ',
       'Mutation Description CDS', 'Mutation Description AA',
       'ONTOLOGY_MUTATION_CODE', 'GENOMIC_MUTATION_ID',
       'Mutation genome position GRCh37', 'Mutation genome position GRCh38',
       'COSMIC_SAMPLE_TESTED', 'COSMIC_SAMPLE_MUTATED', 'DISEASE',
       'WGS_DISEASE', 'EXAC_AF', 'EXAC_AFR_AF', 'EXAC_AMR_AF', 'EXAC_ADJ_AF',
       'EXAC_EAS_AF', 'EXAC_FIN_AF', 'EXAC_NFE_AF', 'EXAC_SAS_AF',
       'GNOMAD_EXOMES_AF', 'GNOMAD_EXOMES_AFR_AF', 'GNOMAD_EXOMES_AMR_AF',
       'GNOMAD_EXOMES_ASJ_AF', 'GNOMAD_EXOMES_EAS_AF', 'GNOMAD_EXOMES_FIN_AF',
       'GNOMAD_EXOMES_NFE_AF', 'GNOMAD_EXOMES_SAS_AF', 'GNOMAD_GENOMES_AF',
       'GNOMAD_GENOMES_AFR_AF', 'GNOMAD_GENOMES_AMI_AF',
       'GNOMAD_GENOMES_AMR_AF', 'GNOMAD_GENOMES_ASJ_AF',
       'GNOMAD_GENOMES_EAS_AF', 'GNOMAD_GENOMES_FIN_AF',
       'GNOMAD_GENOMES_MID_AF', 'GNOMAD_GENOMES_NFE_AF',
       'GNOMAD_GENOMES_SAS_AF', 'CLINVAR_CLNSIG', 'CLINVAR_TRAIT', 'GERP++_RS',
       'MIN_SIFT_SCORE', 'MIN_SIFT_PRED', 'DNDS_DISEASE_QVAL_SIG',
       'MUTATION_SIGNIFICANCE_TIER']

AA = {'ALA':'A',
'ARG': 'R',
'ASN':'N',
'ASP':'D',
'ASX':'B',
'CYS':'C',
'GLU':'E',
'GLN':'Q',
'GLX':'Z',
'GLY':'G',
'HIS':'H',
'ILE':'I',
'LEU':'L',
'LYS':'K',
'MET':'M',
'PHE':'F',
'PRO':'P',
'SER':'S',
'THR':'T',
'TRP':'W',
'TYR':'Y', 
'VAL':'V'
}


#samp = pd.read_csv('/home/jovyan/session_data/mounted-data/LP3000429-DNA_G03.vcf', comment='#', sep='\t', header=None)
samp = pd.read_csv(annotation_vcf_path, comment='#', sep='\t', header=None)
samp = samp.rename(columns={0: 'chr', 1:'pos', 2:'ID', 3:'REF', 4:'ALT', 5:'QUAL', 6:'FILTER', 7:'INFO', 8:'FORMAT',9:'NORMAL', 10:'TUMOR'})
#CHROM POSIDREF ALT QUAL FILTER INFO FORMAT
set(samp['chr'])

unique_chr = pd.DataFrame(set(samp['chr']), columns=['unique_chr'])
unique_chr.to_csv(sample + '_unique_chr.csv', index=False)
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
samp = samp[samp['chr'].isin(chroms)]
samp.to_csv(sample+ '_unfiltered_variants.csv', index=False)
samp_filtered = samp[samp['FILTER']=='PASS']
samp_filtered.to_csv(sample+ '_filtered_variants.csv', index=False)
#samp
##keep rows which contain the relevant terms anywhere in the info column
#to_keep = list(samp['INFO'].str.contains('|'.join(relevant_terms)))
#samp['relevant_types']=to_keep
#sampcsqt_type = samp.loc[samp['relevant_types']==True]
samp_tert = samp[(samp['pos'].isin([1295228, 1295250,1295229,1295242,1295243])) & (samp['chr']=='chr5')]
samp_tert.to_csv(sample + '_tert_promoter_mutations_from_unfiltered_mutation_file.csv', index=False)
