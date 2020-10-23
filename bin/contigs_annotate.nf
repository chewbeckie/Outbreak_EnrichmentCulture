#!/usr/bin/env nextflow

/*
* created by Johanna Wong Oct 2020
*/

/*Parameters*/
params.name //user-defined contigs name (for output)
params.contigs //path to fasta file of contigs
params.krakendb //path to MiniKraken2 database (must be absolute path)
params.plasmidlist //index csv to path of the plasmid genomes database fasta files (must be absolute path)
params.chromosomelist //index csv to path of the chromosome genomes database fasta files (must be absolute path)
params.contigs_as_ref //use this variable if want to use contigs as reference instead of query for mash search

/*Channels*/
contigs = file(params.contigs)

Channel.fromPath(params.plasmidlist)
        .splitText()
        .map { it.replaceFirst(/\n/,'') }
        .into{ plasmids_ch}

Channel.fromPath(params.chromosomelist)
        .splitText()
        .map { it.replaceFirst(/\n/,'') }
        .into{chromosomes_ch}

/*Processes*/

process kraken {
    tag "${contigs}"
    publishDir "./${params.name}_results/kraken", mode: 'copy'
    conda 'bioconda::kraken2'

    output:
        file("*")
        file("*.krakenout") into krakenout_ch
    
    script:
    """
    kraken2 --db $params.krakendb $contigs\
    --unclassified-out ${params.name}_unclassified.tsv --classified-out ${params.name}_classified.tsv \
    --memory-mapping --use-names \
    --report ${params.name}_report.tsv --output ${params.name}.krakenout
    """
    }
    
process mash_plasmids {
    tag "${contigs}"
    publishDir "${params.name}_results", mode: 'copy'
    conda 'bioconda::mash'

    input:
        path(plasmids) from plasmids_ch.toList()

    output:
        file("*.msh")
        file("*.tsv")

    script:
    if (params.contigs_as_ref){
        """
        mash sketch -i $contigs -o ${params.name}_ref.msh
        mash dist ${params.name}_ref.msh $plasmids > ${params.name}_pls_mash_rev.tsv
        """
    } else {
        """
        mash sketch -i $plasmids -o pls_ref.msh
        mash dist pls_ref.msh $contigs > ${params.name}_pls_mash.tsv
        """
    }
}

process mash_chromosomes {
    tag "${contigs}"
    publishDir "${params.name}_results", mode: 'copy'
    conda 'bioconda::mash'

    input:
        path(chromosomes) from chromosomes_ch.toList()

    output:
        file("*.msh")
        file("*.tsv")

    script:
    if (params.contigs_as_ref){
        """
        mash sketch -i $contigs -o ${params.name}_ref.msh
        mash dist ${params.name}_ref.msh $chromosomes > ${params.name}_chr_mash_rev.tsv
        """
    } else {
        """
        mash sketch -i $chromosomes -o chr_ref.msh
        mash dist chr_ref.msh $contigs > ${params.name}_chr_mash.tsv
        """
    }
}

process abricate {
    tag "${contigs}"
    publishDir "${params.name}_results", mode: 'copy'
    conda 'bioconda::abricate'

    output:
        file("*.tsv")
    
    script:
    """
    abricate --db resfinder $contigs > ${params.name}_resfinder_amr.tsv
    abricate --db plasmidfinder $contigs > ${params.name}_plsfinder_pmlst.tsv
    """
}

/***need to add step to summarize the results together
/*process result_collage{
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