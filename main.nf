#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{sample, mutations, regions -> [sample, file(mutations), file(regions)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  extract_coding_mutations{
    //maxForks 900
    //errorStrategy 'finish'
    //maxRetries 0
    tag"$sample"
    publishDir "${params.outdir}/$sample", mode: 'copy'
    
    input:
    set val(sample), file(mutations),file(regions) from ch_input

    output:
    file '*_coding_mutations.csv'
    
    script:
    """
    extract_coding_mutations.py -sample '$sample' -mutations '$mutations' -regions '$regions'
    """ 
}
