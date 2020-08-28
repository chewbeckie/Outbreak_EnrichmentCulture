#!/usr/bin/env nextflow

/*
 * created by Johanna Wong 2 April 2020
 */

/*
Parameter
/*
//run table files with file name and the file to read1 and read2*/
params.index = "/shared/homes/149306/megahit_testing/index_EW_08042020.csv"
//setting output folder locations
params.outdir = "/shared/homes/149306/megahit_testing/output"

//setting for read trimming
params.mean_quality = 15
params.trimming_quality = 20

// Channels

Channel.fromPath(params.index)
        .splitCsv(header:true, sep: ',')
        .map{row-> tuple(row.sampleId, file(row.read1), file(row.read2))}
        .into { samples_ch; trim_ch }

// Step 1 - FASTQC
    process fastqc {
     queue 'smallq'
     memory '5 GB'
     cpus  4
     executor 'pbspro'

        conda 'bioconda::fastqc'
        tag "$sampleId"
        publishDir "${params.outdir}/", mode: 'copy',
            saveAs: {filename -> filename.indexOf(".zip") == -1 ? "fastqc/$filename" : null}

        input:
          set sampleId, file(read1), file(read2) from samples_ch

        output:
          file '*_fastqc.{zip,html}' into fastqc_results

        script:
          """
          fastqc -q $read1
          fastqc -q $read2
          """
    }

//Step 2 - read trimming
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
          val qual from params.mean_quality
          val trim_qual from params.trimming_quality

        output:
//          set val(sampleId), file("${sampleId}_merged.fastq.gz") into merged_reads, merged_reads_mapping
          set val(sampleId), file("trimmed_${sampleId}_R1_001.fastq.gz") into trimmed_read1, read1_mapping
          set val(sampleId), file("trimmed_${sampleId}_R2_001.fastq.gz") into trimmed_read2, read2_mapping
          file("fastp.*")

        script:
//      -m --merged_out "${sampleId}_merged.fastq.gz"\
          """
          fastp -q "${qual}" -5 -3 --correction\
          --cut_mean_quality "${trim_qual}"\
          -i "${read1}" -I "${read2}"\
          -o "trimmed_${sampleId}_R1_001.fastq.gz"\
          -O "trimmed_${sampleId}_R2_001.fastq.gz"
          """
}

// Step 3 - Megahit de novo assembly
process megahit_assembly{
  queue 'workq'
  memory '200 GB'
  cpus  24
  executor 'pbspro'

        tag "$sampleId"
        conda 'bioconda::megahit'
        publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
        set val(sampleId), file(tread1) from trimmed_read1
        set val(sampleId), file(tread2) from trimmed_read2

    output:
        set val(sampleId), file("megahit_out/*.contigs.fa") into contigs_mapping

    script:
        """
        megahit -t 24 -1 $tread1 -2 $tread2 \
        --out-prefix $sampleId
        """
}

// Step 4 - Read mapping back to contigs
process mapping{
  queue 'workq'
  memory '20 GB'
  cpus  8
  executor 'pbspro'

        tag "$sampleId"
        conda 'bioconda::bwa; bioconda::samtools'
        publishDir "${params.outdir}/bwa", mode: 'copy'

    input:
        set val(sampleId), file(contig) from contigs_mapping
        set val(sampleId), file(tread1) from read1_mapping
        set val(sampleId), file(tread2) from read2_mapping

    output:
        set file("*.txt"), file("*.bam")  into results

    script:
        """
        bwa index $contig
        bwa mem $contig $tread1 $tread2 | samtools view -S -b - > ${sampleId}_mapped.bam
        samtools flagstat ${sampleId}_mapped.bam > ${sampleId}_flagstat.txt
        """
}
