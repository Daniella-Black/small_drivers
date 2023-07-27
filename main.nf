#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey, somatic_small_variants_annotation_vcf, mane, hgnc, cmc, non_mane_transcripts, cgc-> [tumour_sample_platekey, file(somatic_small_variants_annotation_vcf), file(mane), file(hgnc), file(cmc), file(non_mane_transcripts), file(cgc)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    maxForks 900
    errorStrategy 'ignore'
    maxRetries 3
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}/$tumour_sample_platekey", mode: 'copy'
    
    input:
    set val(tumour_sample_platekey), file(somatic_small_variants_annotation_vcf), file(mane), file(hgnc), file(cmc), file(non_mane_transcripts), file(cgc) from ch_input

    output:
    file '*_coding_mutations.csv'
    file '*_coding_mutations_no_mane_transcript_maps.csv'
    file '*_coding_mutations_more_than_1_mane_transcript_maps.csv'
    file '*_sampcsqt_type.csv'
    file '*_cosmic_annotation_using_pos_no_indels.csv'
    file '*_one_cgc_mane_tran_associated_with_relevant_terms.csv'
    file '*_over_one_cgc_mane_tran_associated_with_relevant_term.csv'
    file '*_tert_promoter_mutations.csv'
    //file '*_pre_cosmic_for_splice_investigation.csv'
    
    script:
    """
    coding_mutations_nf.py -sample '$tumour_sample_platekey' -annotation_vcf_path '$somatic_small_variants_annotation_vcf' -mane '$mane' -hgnc '$hgnc' -cmc '$cmc' -non_mane_transcripts '$non_mane_transcripts' -cgc '$cgc'
    """ 
}
