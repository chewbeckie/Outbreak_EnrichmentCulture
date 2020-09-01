#!/usr/bin/env nextflow

/*
 * created by Johanna Wong 31 Jul 2020
 */

/*
Parameter
/*
//run table files with sample name and paths to read1, read2 */
params.shortreads //= "$params.InputDir/shortreads_SH.csv" or "$params.InputDir/shortreads_WG.csv"
//sample name and path to long read files
params.longreads //= "$params.InputDir/longreads_SH.csv" or "$params.InputDir/longreads_WG.csv"

//name of contig to be made
params.contigname = //"SH" or "WG"

//setting for metaflye
params.gs = "100m"

// Channels

Channel.fromPath(params.shortreads)
        .splitCsv(header:true, sep: ',')
        .map{row-> tuple(row.sampleId, file(row.read1), file(row.read2))}
        .set { trim_ch }

Channel.fromPath(params.longreads)
        .splitCsv(header:true, sep: ',')
        .map{row-> tuple(row.sampleId, file(row.longread))}
        .set { longread_ch }


//Step 1 - read trimming
process fastp {
  tag "$sampleId"
  publishDir "${params.OutputDir}/trimmed", mode: 'copy'

  input:
    set sampleId, file(read1), file(read2) from trim_ch

  output:
    file("trimmed_${sampleId}_R1.fastq.gz") into tread1_ch
    file("trimmed_${sampleId}_R2.fastq.gz") into tread2_ch
    file("fastp.*") into fastpresult

  script:
    """
    fastp -i "${read1}" -I "${read2}"\
      -o "trimmed_${sampleId}_R1.fastq.gz"\
      -O "trimmed_${sampleId}_R2.fastq.gz"
    """
}

// Step 2 - merging read1 and read2 together

process merge{
  publishDir "${params.OutputDir}/trimmed", mode: 'copy'

  input:
    path '*_R1.fastq.gz' from tread1_ch.toList()
    path '*_R2.fastq.gz' from tread2_ch.toList()
  
  output:
    set file("*_R1.fastq.gz"), file("*_R2.fastq.gz") into mergedreads_ch

  script:
    """
    cat *_R1.fastq.gz > merged_${params.contigname}_R1.fastq.gz
    cat *_R2.fastq.gz > merged_${params.contigname}_R2.fastq.gz
    """

}

// Step 2 - Metaflye assembly
process metaflye_assembly{
  tag "$sampleId"
  publishDir "${params.OutputDir}/assembly", mode: 'copy'

  input:
    set val(sampleId), file(longread) from longread_ch

  output:
    set val(sampleId), file("flye_${sampleId}/*") into assembly_result
    set val(sampleId), file("${sampleId}_assembly.gfa") into graph_ch
    
  script:
    """
    flye --nano-raw $longread \
      --meta --genome-size $params.gs \
      -t 16 -o flye_$sampleId
    cp flye_${sampleId}/assembly_graph.gfa ${sampleId}_assembly.gfa
    """
}

// Step 4 - Polish with ntEdit
process ntedit_polish{  
  queue 'workq'
  memory '40 GB'
  cpus  16
  executor 'pbspro'

        tag "$sampleId"
        conda 'bioconda::ntedit'
        publishDir "${params.OutputDir}/ntedit_polish", mode: 'copy'

    input:
        set file(merged_read1), file(merged_read2) from mergedreads_ch
        set val(sampleId), file(graph) from graph_ch

    output:
        file("*.fa")
        file("*.tsv")
        file("*.vcf")
        file("*.bf")
        file("*.gfa")

    script:
      """
      flye_gfa_to_fasta_for_polishing.py $graph ${sampleId}_assembly.fa
      nthits -c 1 --outbloom -p solidBF_${sampleId} -b 36 -k 40 -t 16 $merged_read1 $merged_read2
      ntedit -m 1 -f ${sampleId}_assembly.fa -r solidBF_${sampleId}_k40.bf -b ${sampleId}_assembly
      flye_polished_fa_to_gfa.py $graph ${sampleId}_assembly_edited.fa polished_${sampleId}.gfa
      """
}
