#!/usr/bin/env nextflow

/*
* created by Johanna Wong Oct 2020
*/

/*Parameters*/
params.genome_list //accession number of genomes to download as database
params.contigs //path to fasta file of contigs for mash search
//params.info //information of the genomes downloaded

/*Channels*/
Channel.fromPath(params.contigs)
        .set{ contigs_ch }

Channel.fromPath(params.genome_list)
        .splitText()
        .map { it.replaceFirst(/\n/,'') }
        .into{ chr_ch ; count_ch}

genome_count = count_ch.count()
//info = file(params.info)

/*Processes*/

process get_fasta {
    tag "$chr_id"
    publishDir "genome_fasta", mode: 'copy'
    conda 'bioconda::entrez-direct'

    input:
        val(chr_id) from chr_ch
    
    output:
        file("${chr_id}.fa") into genomes_ch
        file("log.txt") into done_ch
    
    script:
    """
    efetch -db sequences -format fasta -id $chr_id > ${chr_id}.fa
    echo "download completed" > log.txt
    """
}

//pace mash after all genome download completed
if(done_ch.count() == genome_count ){ready_flag = true} else {ready_flag = false}

process mash_sketch {
    publishDir "sketches", mode: 'copy'
    conda 'bioconda::mash'

    input:
        file(contigs) from contigs_ch
    
    when:
        ready_flag = true
    
    output:
        file("*.msh") into sketch_ch
  
    script:
    """
    mash sketch -i $contigs
    """
}

process mash_dist {
    publishDir "mash_dist", mode: 'copy'
    conda 'bioconda::mash'

    input:
        path(ref) from sketch_ch
        file(query) from genomes_ch
    
    output:
        file("*.tsv") into mash_result

    script:
    """
    mash dist -v 0.1 -d 0.1 $ref $query > ${query.baseName}_mash.tsv
    """
}

/*process result_collate{
    publishDir "mash_dist", mode:'copy'

    input:
        file(mash_out) from mash_result
    
    output:
        file("*.tsv")

    script:
    """
    chr_mash_result_edit.R $mash_out $info signf_mash_chr.tsv
    """

}*/