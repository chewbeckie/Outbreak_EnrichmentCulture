#!/usr/bin/env nextflow

/*
* created by Johanna Wong Oct 2020
*/

/*Parameters*/
params.genome_list //accession number of genomes to download as database
params.query //path to fasta file of contigs for mash search

/*Channels*/
Channel.fromPath(params.query)
        .set{ query_ch}

Channel.fromPath(params.genome_list)
        .splitText()
        .map { it.replaceFirst(/\n/,'') }
        .into{ chr_ch ; count_ch}

genome_count = count_ch.count()


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
    publishDir "genome_sketches", mode: 'copy'
    conda 'bioconda::mash'

    input:
        file(genome) from genomes_ch.toList()
    
    when:
        ready_flag = true
    
    output:
        file("*.msh") into sketch_ch
  
    script:
    """
    mash sketch $genome -o genome_collection.msh
    """
}

process mash_dist {
    publishDir "mash_dist", mode: 'copy'
    conda 'bioconda::mash'

    input:
        path(ref) from sketch_ch
        file(query) from query_ch
    
    output:
        file("*.tsv")

    script:
    """
    mash dist -i $ref $query > ${query}_mash.tsv
    """
}
