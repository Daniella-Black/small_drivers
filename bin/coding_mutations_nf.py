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
                       help='path to cgc table)
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
mane[['transcript_ID', 'to_del2']] = mane['name'].str.split('.', 1, expand=True)
mane[['gene_ID', 'to_del2']] = mane['geneName'].str.split('.', 1, expand=True)
mane = mane.rename(columns={'#chrom': 'chr', 'chromStart': 'start', 'chromEnd': 'end', 'geneName2': 'gene_name' })
mane['chr'] = mane['chr'].str.replace('chr', '')
mane_full = mane                       
mane = mane[['chr', 'start', 'end', 'transcript_ID','gene_ID', 'gene_name']]

non_mane_transcripts = pd.read_csv(non_mane_transcripts)
mane = pd.concat([mane, non_mane_transcripts])
mane = mane.reset_index(drop = True)

 

#check how many mane canonical transcripts are listed in the info field. One is good. Sometimes there is more than one. Sometimes there are none. Sometimes there are more than one check_transcripts = list()
#provide a table with an INFO column from which you want to pull out one of a specific set of transcripts of interest e.g. sampcsqt_type
#provide a table with column transcript ID which contains the set of transcripts you are interested in pulling out e.g. mane
##funciton will add a column to the table with the INFO colum in showing which transcripts the list have been identified in the INFO column
def pull_out_transcripts_from_info(table, transcripts_of_interest):
  check_transcripts = list()
  for gene in range(len(table)):
      total_mane_tran_in_row=0
      check_transcripts_tmp =[]
      for transcript in list(transcripts_of_interest['transcript_ID']):
          if transcript in table['INFO'][gene]:
              check_transcripts_tmp.append(transcript)
              total_mane_tran_in_row = total_mane_tran_in_row+1
      if total_mane_tran_in_row == 0:
          check_transcripts.append(['NoManeTran'])
      else:
          check_transcripts.append(check_transcripts_tmp) 

  table['mane_tran'] = check_transcripts
  return table



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

cmc= pd.read_csv(cmc, sep=',', header = None, names=cmc_colnames)

cmc= cmc[['GENE_NAME','ACCESSION_NUMBER','ONC_TSG','CGC_TIER','MUTATION_URL','LEGACY_MUTATION_ID','Mutation CDS','Mutation Description CDS','Mutation Description AA','Mutation AA', 'AA_MUT_START',
       'AA_MUT_STOP', 'SHARED_AA', 'GENOMIC_WT_ALLELE_SEQ',
       'GENOMIC_MUT_ALLELE_SEQ', 'AA_WT_ALLELE_SEQ', 'AA_MUT_ALLELE_SEQ','ONTOLOGY_MUTATION_CODE',
'GENOMIC_MUTATION_ID','Mutation genome position GRCh38','COSMIC_SAMPLE_TESTED','COSMIC_SAMPLE_MUTATED']]

hgnc = pd.read_csv(hgnc,sep='\t')



#samp = pd.read_csv('/home/jovyan/session_data/mounted-data/LP3000429-DNA_G03.vcf', comment='#', sep='\t', header=None)
samp = pd.read_csv(annotation_vcf_path, comment='#', sep='\t', header=None)
samp = samp.rename(columns={0: 'chr', 1:'pos', 2:'ID', 3:'REF', 4:'ALT', 5:'QUAL', 6:'FILTER', 7:'INFO', 8:'FORMAT',9:'NORMAL', 10:'TUMOR'})
#CHROM POSIDREF ALT QUAL FILTER INFO FORMAT
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
samp = samp[samp['chr'].isin(chroms)]
samp = samp[samp['FILTER']=='PASS']
#samp
##keep rows which contain the relevant terms anywhere in the info column
to_keep = list(samp['INFO'].str.contains('|'.join(relevant_terms)))
samp['relevant_types']=to_keep
sampcsqt_type = samp.loc[samp['relevant_types']==True]
samp_tert = samp[(samp['pos'].isin([1295228, 1295250,1295229,1295242,1295243])) & (samp['chr']=='chr5')]
samp_tert.to_csv(sample + '_tert_promoter_mutations.csv', index=False)
##join together variants which have a mutation in the important TERT promoter sites and those which are associated with relevant terms
#sampcsqt_type = pd.concat([sampcsqt_type, samp_tert])
sampcsqt_type.index = pd.RangeIndex(len(sampcsqt_type.index))


#check how many mane canonical transcripts are listed in the info field. One is good. Sometimes there is more than one. Sometimes there are none. Sometimes there are more than one check_transcripts = list()
pull_out_transcripts_from_info(sampcsqt_type, mane)



##output these tables to investigate strange cases
sampcsqt_type_over_1 = sampcsqt_type.loc[sampcsqt_type['mane_tran'].map(len)>1]
sampcsqt_type_over_1.index = pd.RangeIndex(len(sampcsqt_type_over_1.index))
sampcsqt_type_over_1.to_csv(sample + '_coding_mutations_more_than_1_mane_transcript_maps.csv')
sampcsqt_type_nocantran =sampcsqt_type.loc[sampcsqt_type['mane_tran'].map(len)==1]
#filter out rows where there are no canonical transcripts in the info column 
sampcsqt_type_nocantran['mane_tran']=flatten(sampcsqt_type_nocantran['mane_tran'])
sampcsqt_type_nocantran = sampcsqt_type_nocantran.loc[sampcsqt_type_nocantran['mane_tran'] == 'NoManeTran']
sampcsqt_type_nocantran.to_csv(sample + '_coding_mutations_no_mane_transcript_maps.csv')
###

#filter out rows where there is more than one canonical transcript in the info column
sampcsqt_type = sampcsqt_type.loc[sampcsqt_type['mane_tran'].map(len)==1]
sampcsqt_type['mane_tran'] = flatten(sampcsqt_type['mane_tran'])
sampcsqt_type = sampcsqt_type.loc[sampcsqt_type['mane_tran'] != 'NoManeTran']
sampcsqt_type.index = pd.RangeIndex(len(sampcsqt_type.index))
sampcsqt_type.to_csv(sample + '_sampcsqt_type.csv')

if len(sampcsqt_type_over_1.index) >0:
    variant_info = list()
    variant_info_transcript = list()
    VAF=list()
    for row in range(len(sampcsqt_type_over_1.index)): 
        variant_info_temp = list()
        enst_is_coding = list()
        enst_is_coding_transcript = list()
        if sampcsqt_type_over_1['INFO'][row].find('VAF')!= -1:
            VAF.append(sampcsqt_type_over_1['INFO'][row].rsplit('VAF=', maxsplit=1)[1].rsplit(';')[0])
        else:
            VAF.append('No Value')
        for gene in sampcsqt_type_over_1['mane_tran'][row]:
            if len(re.findall (gene, sampcsqt_type_over_1['INFO'][row]))==1: 
                #split at the ENST to next,|;
                outtemp = sampcsqt_type_over_1['INFO'][row].rsplit(gene, maxsplit=1)[1]
                outtemp = re.split(',|;', outtemp)[0]
                variant_info_temp.append(outtemp + '$' + gene)
            if len(re.findall (gene, sampcsqt_type_over_1['INFO'][row])) ==2:
                #split at both ENST to next,|;
                outtemp = sampcsqt_type_over_1['INFO'][row].rsplit(gene, maxsplit=2)[2]
                outtemp = re.split(',|;', outtemp)[0]
                variant_info_temp.append(outtemp+ '$' + gene)
        for info in variant_info_temp:
            if any(substring in info for substring in relevant_terms):
                enst_is_coding.append(info)
                enst_is_coding_transcript.append(gene)
        variant_info.append(list(enst_is_coding))
        variant_info_transcript.append(list(enst_is_coding_transcript))
    sampcsqt_type_over_1['relevant_term_associated_trans'] = variant_info_transcript
    for element in variant_info:
        if len(element) >1:
            if element[0] == element[1]:
                element.remove(element[1])
    sampcsqt_type_over_1['variant_info'] =variant_info
    sampcsqt_type_over_1['VAF'] =VAF
    sampcsqt_type_over_1_multi =sampcsqt_type_over_1.loc[sampcsqt_type_over_1['variant_info'].map(len)>1].reset_index(drop=True)
    sampcsqt_type_over_1=sampcsqt_type_over_1.loc[sampcsqt_type_over_1['variant_info'].map(len)==1].reset_index(drop=True)
    if len(sampcsqt_type_over_1.index) > 0:
        sampcsqt_type_over_1['variant_info'] = flatten(sampcsqt_type_over_1['variant_info'])
        sampcsqt_type_over_1[['variant_info', 'mane_tran']] = sampcsqt_type_over_1['variant_info'].str.split('$', 1, expand=True)
    else:
        sampcsqt_type_over_1 = None
                       
##################################################################################################################################
#################dealing with cases where there are more than one mane transcript associated with relevant term in the info column
##################################################################################################################################                       
##find out how many cgc genes are in the list of mane transcripts associated with relevant terms for the variants with multiple mane transcripts associated with relevant terms in the info column
mane_cgc = mane[mane['gene_ID'].isin(list(cgc['ENSG']))].reset_index(drop=True)
at_least_one_tran_cgc = list()
for row in range(len(sampcsqt_type_over_1_multi)):
    at_least_one_tran_cgc = 0
    for gene in range(len(mane_cgc)):
        if mane_cgc['transcript_ID'][gene] in sampcsqt_type_over_1_multi['relevant_term_associated_trans'][row]:
            at_least_one_tran_cgc = at_least_one_tran_cgc +1
    if at_least_one_tran_cgc == 1:
        at_least_one_tran_cgc.apend('1_cgc')
    elif at_least_one_tran_cgc == 2:
        at_least_one_tran_cgc.apend('2_cgc')
    elif at_least_one_tran_cgc == 3:
        at_least_one_tran_cgc.apend('3_cgc')
    elif at_least_one_tran_cgc > 3:
        at_least_one_tran_cgc.apend('over_3_cgc')
    else:
        at_least_one_tran_cgc.apend('0_cgc')
sampcsqt_type_over_1_multi['cgc_tran_count'] = at_least_one_tran_cgc

##get strand information for all the mane transcripts associated with relevant terms overlapping variant
both_on_strand = list()
big_strand = list()
for row in range(len(sampcsqt_type_over_1_multi.index)):
    strand = list()
    for enst in sampcsqt_type_over_1_multi['relevant_term_associated_trans'][row]:
        mane_full_tmp = mane_full[mane_full['transcript_ID']==enst]
        mane_full_tmp = mane_full_tmp.reset_index(drop=True)
        if len(mane_full_tmp.index)>0:
            strand.append(mane_full_tmp['strand'][0])
    big_strand.append(strand)
 
##if one mane tran associated with relevant term is a cgc and the others are not. Then extract the information for the cgc gene and append this to the mane table
sampcsqt_type_over_1_cgc1 = sampcsqt_type_over_1_multi[sampcsqt_type_over_1_multi['cgc_tran_count']=='1_cgc'] 
##if more than one mane tran associated with relevant term is cgc and the others are not. Output for now.
sampcsqt_type_over_1_cgcover1 = sampcsqt_type_over_1_multi[sampcsqt_type_over_1_multi['cgc_tran_count'].isin('2_cgc', '3_cgc', 'over_3_cgc')]                        
                       


#pull out the variant info - split into the two cases os whether the ENST is written once or twice in the info column SomaticFisherPhred = list()
if len(sampcsqt_type.index) >0:
  VAF= list()
  phyloP = list()
  variant_info = list()
  for gene in range(len(sampcsqt_type.index)):
      if sampcsqt_type['INFO'][gene].find('VAF') != -1:
          VAF.append(sampcsqt_type['INFO'][gene].rsplit('VAF=', maxsplit=1)[1].rsplit(';')[0])
      else:
          VAF.append('No Value')
      if len(re.findall(sampcsqt_type['mane_tran'][gene], sampcsqt_type['INFO'][gene])) ==1:
          #split at the ENST to next,|;
          outtemp = sampcsqt_type['INFO'][gene].rsplit(sampcsqt_type['mane_tran'][gene], maxsplit=1)[1]
          outtemp = re.split(',|;', outtemp)[0]
          variant_info.append(outtemp)
      if len(re.findall(sampcsqt_type['mane_tran'][gene], sampcsqt_type['INFO'][gene])) ==2:
          #split at both ENST to next,];
          outtemp = sampcsqt_type['INFO'][gene].rsplit(sampcsqt_type['mane_tran'][gene], maxsplit=2)[2]
          outtemp = re.split(',|;', outtemp)[0]
          variant_info.append(outtemp)


  sampcsqt_type['variant_info'] = variant_info
  sampcsqt_type['VAF']=VAF
if len(sampcsqt_type.index)>0:
  if sampcsqt_type_over_1 is not None and len(sampcsqt_type_over_1.index)>0:
          sampcsqt_type_full = pd.concat([sampcsqt_type,sampcsqt_type_over_1])
          sampcsqt_type_full = sampcsqt_type_full[['chr', 'pos', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'mane_tran', 'variant_info', 'VAF']] 
          sampcsqt_type_full.to_csv(sample + '_coding_mutations.csv', index=False)
          coding = sampcsqt_type_full
  else:
    sampcsqt_type = sampcsqt_type[['chr', 'pos', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'mane_tran', 'variant_info', 'VAF']] 
    sampcsqt_type.to_csv(sample + '_coding_mutations.csv', index=False)
    coding = sampcsqt_type
elif len(sampcsqt_type.index)<1 and len(sampcsqt_type_over_1.index)>0:
    sampcsqt_type_over_1 = sampcsqt_type_over_1[['chr', 'pos', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'mane_tran', 'variant_info', 'VAF']] 
    sampcsqt_type_over_1.to_csv(sample + '_coding_mutations.csv', index=False)
    coding = sampcsqt_type_over_1
else:
  df = pd.DataFrame(columns=['chr', 'pos', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'mane_tran', 'variant_info', 'VAF'])
  df.to_csv(sample + '_coding_mutations.csv', index=False)
  coding = df
  

sampcsqt_type = None
sampcsqt_type_over_1 = None
samp = None
sampcsqt_type_full = None
variant_info = None

if len(coding.index) > 0:
  #############
  ##with hgnc redundant names
  #####
  #################
  #prepare the hgnc database to get all the gene names/previoyus names and aliasses in a list
  #hgnc = pd.read_csv('/home/jovyan/session_data/MutationTimeR_Gel/hgnc_table_with_locus_type_and_ensembl_gene_id.txt', sep='\t') 
  hgnc = hgnc[hgnc['Locus type'] == 'gene with protein product']
  hgnc['all_names'] = hgnc['Approved symbol'].astype('str') + ',' + hgnc['Previous symbols'].astype('str')+','+ hgnc['Alias symbols'].astype('str') 
  hgnc['all_names'] = hgnc['all_names'].str.replace(',nan', '')
  hgnc = hgnc.reset_index(drop=True)
  for row in range(len(hgnc.index)):
      hgnc['all_names'][row] = hgnc['all_names'][row].split(',')
  #dictionary of redundant gene names and ENSG
  hgnc_d = dict(zip(hgnc['Ensembl gene ID'], hgnc['all_names']))
  #make a dictionary of the ENSG names for each GENE_NAME in cosmic by searching against the redundant gene_names dictionary 
  cosmic_dict = {}
  for cosmic_gene in set(cmc['GENE_NAME']):
      for ENSG in hgnc_d:
          if cosmic_gene in hgnc_d[ENSG]:
              cosmic_dict[cosmic_gene] = ENSG
  ##convert to a dataframe and merge with the full cosmic database so that the ENSG names are available for every cosmic entry 
  cosmic_ENSG_GENE_NAME = pd.DataFrame.from_dict(cosmic_dict, orient='index') 
  cosmic_ENSG_GENE_NAME = cosmic_ENSG_GENE_NAME.reset_index()
  cosmic_ENSG_GENE_NAME = cosmic_ENSG_GENE_NAME.rename(columns={0:'ENSG', 'index': 'GENE_NAME'}) #cosmic_ENSG_GENE_NAME
  cosmic_ENSG = cmc.merge(cosmic_ENSG_GENE_NAME, on="GENE_NAME")


  coding[['empty', 'cds', 'aa', 'effect', 'remaining_variant_info']] = coding['variant_info'].str.split('|', 4, expand=True)
  coding = coding[coding['effect'].isin(relevant_terms)]
  coding_nonspliceregion = coding[coding['effect'] != 'splice_region_variant']
  coding_spliceregion = coding[coding['effect'] == 'splice_region_variant']
  coding_spliceregion = coding_spliceregion[coding_spliceregion['INFO'].str.contains('\+5(?=[A-Z])')] 
  coding = pd.concat([coding_nonspliceregion,coding_spliceregion])
  coding.index = pd.RangeIndex(len(coding.index))
  for amino in AA.keys():
      for row in range(len(coding.index)):
          if amino in coding['aa'][row]:
              coding['aa'][row] = coding['aa'][row].replace(amino, AA[amino])

  ##split out the subs and indels because are dealing with these differently for cosmic. 
  coding['aa'] = coding['aa'].str.replace(r'(', '', regex=True)
  coding['aa'] = coding['aa'].str.replace(r')', '', regex=True)
  #coding_indels = coding.loc[(coding['REF'].str.len() > 1) | (coding['ALT'].str.len() >1)]
  #coding_indels.to_csv('/home/jovyan/session_data/small_drivers/Head_Neck/mane/including_over_1/LP3001647-DNA_A05/LP3001647-DNA_A05_coding_mutations_indels_one_letter_code.csv')
  #coding = coding.loc[(coding['REF'].str.len() == 1) & (coding['ALT'].str.len() == 1) ]
  coding.index = pd.RangeIndex(len(coding.index))

  for row in range(len(coding.index)):
      if 'STOP' in coding['aa'][row]:
          coding['aa'][row] = coding['aa'][row].replace('STOP', '*')
  for row in range(len(coding.index)):
      if len(coding['REF'][row]) == 1 and len(coding['ALT'][row]) == 1:
          for amino in AA.values():
              if coding['aa'][row].count(amino) ==2:
                  coding['aa'][row] = nth_repl(coding['aa'][row], amino,"=",2) #print('yes')
  ##make the pos column
  pos_for_cosmic_comp = list()
  for row in range(len(coding.index)):
      if len(coding['REF'][row]) == 1 and len(coding['ALT'][row])==1:
          if coding['chr'][row] == 'chrX':
            pos_for_cosmic_comp.append('23:' + str(coding['pos'][row]) +'-' + str(coding['pos'][row]))
          elif coding['chr'][row] == 'chrY':
            pos_for_cosmic_comp.append('24:' + str(coding['pos'][row]) +'-' + str(coding['pos'][row]))
          else:
            pos_for_cosmic_comp.append(coding['chr'][row] + ':' + str(coding['pos'][row]) +'-' + str(coding['pos'][row]))
      elif len(coding['REF'][row]) == 1 and len(coding['ALT'][row]) > 1:
          pos_for_cosmic_comp.append('indel')
          #pos_for_cosmic_comp.append(coding['chr'][row] + ':' + str(coding['pos'][row]) +'-' + str(coding['pos'][row] + 1))
      elif len(coding['REF'][row]) > 1 and len(coding['ALT'][row]) == 1:
          pos_for_cosmic_comp.append('indel')
          #pos_for_cosmic_comp.append(coding['chr'][row] + ':' + str(coding['pos'][row]) +'-' + str(coding['pos'][row] + len(coding['REF'][row])- 1))
      else:
          pos_for_cosmic_comp.append('strange_indel')

  coding['pos_for_cosmic_comp'] = pos_for_cosmic_comp
  coding['pos_for_cosmic_comp'] = coding['pos_for_cosmic_comp'].str.replace('chr', '')
  
  #coding.to_csv(sample+ '_pre_cosmic_for_splice_investigation.csv')
  ##gsub the chr
  #combine the man transcript table with the coding mutation table so can convert from ENST to ENSG to match cosmic.
  mane_temp= mane.rename(columns={'chr':'chr_mane', 'start': 'start_mane', 'end': 'end_mane', 'gene_name': 'gene_name_mane', 'transcript_ID':'mane_tran', 'gene_ID': 'gene_ID_mane' })
  mane_temp = mane_temp[mane_temp['mane_tran'].isin(coding['mane_tran'])]
  coding_mane = coding.merge(mane_temp, on='mane_tran')
  #merge cosmic and coding table via ENSG numbers
  coding_mane = coding_mane.rename(columns={'gene_ID_mane': 'ENSG'})
  ##cosmic temp is for creating coding_mane_no_pos where both ENSG and protien must match
  ##for pos can use just the genome position
  cosmic_temp = cosmic_ENSG[cosmic_ENSG['ENSG'].isin(coding_mane['ENSG'])]
  cosmic_pos = cmc[cmc['Mutation genome position GRCh38'] != 'NaN']
  cosmic_temp_no_pos = cosmic_temp[cosmic_temp['Mutation genome position GRCh38'] == 'NaN']

  coding_mane_pos = coding_mane.merge(cosmic_pos, left_on='pos_for_cosmic_comp', right_on='Mutation genome position GRCh38') 
  if len(cosmic_temp_no_pos.index) >0:
      coding_mane_no_pos = coding_mane.merge(cosmic_temp_no_pos, on='ENSG')
      coding_mane_no_pos = coding_mane.loc[(coding_mane['aa'] ==coding_mane['Mutation AA'])]
      #coding_mane_pos = coding_mane_pos.loc[(coding_mane_pos['Mutation genome position GRCh38'] ==coding_mane_pos['pos_for_cosmic_comp'])]
  #select rows which have cosmic info.
  if len(cosmic_temp_no_pos.index) >0:
      coding_cosmic_final = pd.concat[[coding_mane_pos, coding_mane_no_pos]]
  else:
      coding_cosmic_final = coding_mane_pos
  #coding_mane = coding_mane.loc[(coding_mane['aa'] ==coding_mane['Mutation AA'])]
  coding_cosmic_final.to_csv(sample+'_cosmic_annotation_using_pos_no_indels.csv')
else:
  coding_cosmic_final = pd.DataFrame()
  coding_cosmic_final.to_csv(sample+'_cosmic_annotation_using_pos_no_indels.csv')                                                  
                                                    
