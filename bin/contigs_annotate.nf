#!/usr/bin/env nextflow

/*
* created by Johanna Wong Oct 2020
*/

/*Parameters*/
params.name //user-defined contigs name (for output)
params.contigs// ="$baseDir/input_files/SH_assembly_edited.fa" //path to fasta file of contigs
params.krakendb// ="$baseDir/input_files/minikraken2_v1_8GB" //path to MiniKraken2 database (must be absolute path)
params.plspath// ="$baseDir/input_files/plsdb_genomes/*.fa" // path + matching glob pattern of the plasmid genomes database fasta files
params.chrpath// ="$baseDir/input_files/chr_genomes/*.fa" //path + matching glob pattern of the chromosome genomes database fasta files
params.contigs_as_ref //use this variable if want to use contigs as reference instead of query for mash search

/*Channels*/
contigs = file(params.contigs)

//Channel.fromPath(params.plspath)
//        .set{ plasmids_ch }

//Channel.fromPath(params.chrpath)
//        .set{ chromosomes_ch }

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

//Comment:: when dealing with large plasmids/chromosomes database
//it is needed to merge plasmids/chromosomes fasta files into one to avoid argument too long error

process merge_fasta {
    tag "${contigs}"
    publishDir "${params.name}_results", mode: 'copy'

    output:
        file("pls*") into pls_collection
        file("chr*") into chr_collection

    script:
    """
    for f in $params.plspath; do cat \$f >> pls_genome_collection.fa; done
    for f in $params.chrpath; do cat \$f >> chr_genome_collection.fa; done
    """
}


process mash_plasmids {
    tag "${contigs}"
    publishDir "${params.name}_results", mode: 'copy'
    conda 'bioconda::mash'

    input:
        file(pls_db) from pls_collection

    output:
        file("*.msh")
        file("*.tsv")

    script:
    if (params.contigs_as_ref){
        """
        mash sketch -i $contigs -o ${params.name}_ref.msh
        mash dist -i -v 0.1 -d 0.1 ${params.name}_ref.msh $pls_db > ${params.name}_pls_mash_rev.tsv
        """
    } else {
        """
        mash sketch -i $pls_db -o pls_ref.msh
        mash dist -i -v 0.1 -d 0.1 pls_ref.msh $contigs > ${params.name}_pls_mash.tsv
        """
    }
}

process mash_chromosomes {
    tag "${contigs}"
    publishDir "${params.name}_results", mode: 'copy'
    conda 'bioconda::mash'

    input:
        file(chr_db) from chr_collection

    output:
        file("*.msh")
        file("*.tsv")

    script:
    if (params.contigs_as_ref){
        """
        mash sketch -i $contigs -o ${params.name}_ref.msh
        mash dist -i -v 0.1 -d 0.1 ${params.name}_ref.msh $chr_db > ${params.name}_chr_mash_rev.tsv
        """
    } else {
        """
        mash sketch -i $chr_db -o chr_ref.msh
        mash dist -i -v 0.1 -d 0.1 chr_ref.msh $contigs > ${params.name}_chr_mash.tsv
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