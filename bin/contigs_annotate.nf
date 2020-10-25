#!/usr/bin/env nextflow

/*
* created by Johanna Wong Oct 2020
*/

/*Parameters*/
params.name //user-defined contigs name (for output)
params.contigs //path to fasta file of contigs
params.krakendb //path to MiniKraken2 database (must be absolute path)
params.plspath //path + matching glob pattern of the plasmid genomes database fasta files
params.chrpath //path + matching glob pattern of the chromosome genomes database fasta files
params.chr_info //path to chromosome information file e.g. complete_assembly_enterobacteriaceae_chromosome_info.csv
params.pls_info //path to plasmid information file e.g. plsdb.tsv
params.plsmsh //path to plasmid database .msh mash sketch file, if supplied, contigs_as_ref and plspath must be absent
params.chrmsh //path to chromosome database .msh mash sketch file, if supplied, contigs_as_ref and chrpath must be absent
params.contigs_as_ref //use this parameter if want to use contigs as reference instead of query for mash search

/*Channels*/
contigs = file(params.contigs)

/*Processes*/

//Comment:: when dealing with large plasmids/chromosomes database
//it is needed to merge plasmids/chromosomes fasta files into one to avoid argument too long error
process build_plsdb {
    tag "${contigs}"
    publishDir "${params.name}_results", mode: 'copy'
    conda 'bioconda::mash'

    output:
        file("*.fa") into pls_fa
        file("*.msh") into pls_db

    script:
    if(params.plsmsh){
        """
        cp $params.plsmsh pls_ref.msh
        echo "" > empty.fa
        """
    }else{
        """
        for f in $params.plspath; do cat \$f >> pls_genome_collection.fa; done
        mash sketch -i pls_genome_collection.fa -o pls_ref.msh
        """
    }
}

process mash_plasmids {
    tag "${contigs}"
    publishDir "${params.name}_results", mode: 'copy'
    conda 'bioconda::mash'

    input:
        file(msh) from pls_db
        file(fa) from pls_fa
        
    output:
        file("*.tsv")
    
    script:
    if (params.contigs_as_ref){
        """
        mash sketch -i $contigs -o ${params.name}_ref.msh
        mash dist -i -v 0.1 -d 0.1 ${params.name}_ref.msh $fa > ${params.name}_pls_mash_rev.tsv
        """
    } else {
        """
        mash dist -i -v 0.1 -d 0.1 $msh $contigs > ${params.name}_pls_mash.tsv
        """
    }
}

process build_chrdb {
    tag "${contigs}"
    publishDir "${params.name}_results", mode: 'copy'
    conda 'bioconda::mash'

    output:
        file("*.fa") into chr_fa
        file("*.msh") into chr_db

    script:
    if(params.chrmsh){
        """
        cp $params.chrmsh chr_ref.msh
        echo "" > empty.fa
        """
    }else{
        """
        for f in $params.chrpath; do cat \$f >> chr_genome_collection.fa; done
        mash sketch -i chr_genome_collection.fa -o chr_ref.msh
        """
    }
}


process mash_chromosomes {
    tag "${contigs}"
    publishDir "${params.name}_results", mode: 'copy'
    conda 'bioconda::mash'

    input:
        file(msh) from chr_db
        file(fa) from chr_fa
        
    output:
        file("*.tsv")
    
    script:
    if (params.contigs_as_ref){
        """
        mash sketch -i $contigs -o ${params.name}_ref.msh
        mash dist -i -v 0.1 -d 0.1 ${params.name}_ref.msh $fa > ${params.name}_chr_mash_rev.tsv
        """
    } else {
        """
        mash dist -i -v 0.1 -d 0.1 $msh $contigs > ${params.name}_chr_mash.tsv
        """
    }
}

process kraken {
    tag "${contigs}"
    publishDir "./${params.name}_results/kraken", mode: 'copy'
    conda 'bioconda::kraken2'

    output:
        file("*")
        file("*.krakenout") into krakenout_ch //for adding genome_map.py step later
    
    script:
    """
    kraken2 --db $params.krakendb $contigs\
    --unclassified-out ${params.name}_unclassified.tsv --classified-out ${params.name}_classified.tsv \
    --memory-mapping --use-names \
    --report ${params.name}_report.tsv --output ${params.name}.krakenout
    """
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

//summarize the results together (not for contigs_as_ref setting)
process result_collage{
    publishDir "${params.name}_results", mode: 'copy'
    conda 'r::r-tidyverse'

    input:
        file(chr_mash_out) from chr_result
        file(pls_mash_out) from pls_result
    
    output:
        file("*.tsv")

    script:
    if (params.contigs_as_ref){
        """
        chr_mash_result_edit_rev.R $chr_mash_out $params.chr_info ${params.name}_signf_mash_chr_rev.tsv
        plsdb_mash_result_edit_rev.R $pls_mash_out $params.pls_info ${params.name}_signf_mash_pls_rev.tsv
        """
    }else{
        """
        chr_mash_result_edit.R $chr_mash_out $params.chr_info ${params.name}_signf_mash_chr.tsv
        plsdb_mash_result_edit.R $pls_mash_out $params.pls_info ${params.name}_signf_mash_pls.tsv
        """
    }
}