#!/usr/bin/env nextflow

/*
* created by Johanna Wong Oct 2020
*/

/*Parameters*/
params.acc_list //accession number of genomes to download as database
params.out_dir //path to output directory

/*Channels*/
Channel.fromPath(params.acc_list)
        .splitText()
        .map { it.replaceFirst(/\n/,'') }
        .set{ acc_ch }

/*Processes*/

process get_fasta {
    tag "$acc_id"
    publishDir "$params.out_dir", mode: 'copy'
    conda 'bioconda::entrez-direct'

    input:
        val(acc_id) from acc_ch
    
    output:
        file("${acc_id}.fa")
    
    script:
    """
    efetch -db sequences -format fasta -id $acc_id > ${acc_id}.fa
    """
}