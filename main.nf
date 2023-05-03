#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey, somatic_small_variants_annotation_vcf, mane, hgnc, cmc-> [tumour_sample_platekey, file(somatic_small_variants_annotation_vcf), file(mane), file(hgnc), file(cmc)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    errorStrategy 'ignore'
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}/$tumour_sample_platekey", mode: 'copy'
    
    input:
    set val(tumour_sample_platekey), file(somatic_small_variants_annotation_vcf), file(mane), file(hgnc), file(cmc) from ch_input

    output:
    //file '*_all_variant_effects.csv'
    //file '*_splice_variants_of_interest.csv'
    file '*_info_contains_plus5_followed_byA-Z.csv'
    
    script:
    """
    coding_mutations_nf.py -sample '$tumour_sample_platekey' -annotation_vcf_path '$somatic_small_variants_annotation_vcf' -mane '$mane' -hgnc '$hgnc' -cmc '$cmc'
    """ 
}
