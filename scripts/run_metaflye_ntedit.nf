#!/usr/bin/env nextflow

/*
 * created by Johanna Wong 31 Jul 2020
 */

/*
Parameter
/*
//run table files with sample name and paths to read1, read2 */
params.shortreads = "$baseDir/shortreads.csv"
//long read files
params.longreads = "$baseDir/longreads.csv"
//setting output folder locations
params.outdir = "$baseDir/output_files"

//setting for read trimming
params.mean_quality = 15
params.trimming_quality = 20

//setting for metaflye
params.gs = "50m"

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
    queue 'smallq'
    memory '5 GB'
    cpus  4
    executor 'pbspro'

        conda 'bioconda::fastp'
        tag "$sampleId"
        publishDir "${params.outdir}/trimmed", mode: 'copy'

        input:
          set sampleId, file(read1), file(read2) from trim_ch

        output:
          set val(sampleId), file("trimmed_${sampleId}_R1.fastq.gz") into tread1_ch
          set val(sampleId), file("trimmed_${sampleId}_R2.fastq.gz") into tread2_ch
          file("fastp.*") into fastpresult

        script:
          """
          fastp -i "${read1}" -I "${read2}"\
          -o "trimmed_${sampleId}_R1.fastq.gz"\
          -O "trimmed_${sampleId}_R2.fastq.gz"
          """
}

// Step 2 - Metaflye assembly
process metaflye_assembly{
  queue 'workq'
  memory '100 GB'
  cpus  20
  executor 'pbspro'

        tag "$sampleId"
        conda 'bioconda::flye'
        publishDir "${params.outdir}/assembly", mode: 'copy'

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
        publishDir "${params.outdir}/ntedit_polish", mode: 'copy'

    input:
        set val(sampleId), file(tread1) from tread1_ch
        set val(sampleId), file(tread2) from tread2_ch
        set val(sampleId), file(graph) from graph_ch

    output:
        file("*.fa")
        file("*.tsv")
        file("*.vcf")
        file("*.bf")
        file("*.gfa")

/*    shell:
        '''
        awk '/^S/{print ">"$2 "\\n"$3}' !{graph} | fold > !{sampleId}_assembly.fa
        nthits -c 1 --outbloom -p solidBF_!{sampleId} -b 36 -k 40 -t 16 !{tread1} !{tread2}
        ntedit -m 1 -f !{sampleId}_assembly.fa -r solidBF_!{sampleId}_k40.bf -b !{sampleId}
        '''
*/
    script:
      """
      flye_gfa_to_fasta_for_polishing.py $graph ${sampleId}_assembly.fa
      nthits -c 1 --outbloom -p solidBF_${sampleId} -b 36 -k 40 -t 16 $tread1 $tread2
      ntedit -m 1 -f ${sampleId}_assembly.fa -r solidBF_${sampleId}_k40.bf -b ${sampleId}_assembly
      flye_polished_fa_to_gfa.py $graph ${sampleId}_assembly_edited.fa polished_${sampleId}.gfa
      """
}
