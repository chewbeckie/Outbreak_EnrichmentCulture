#!/usr/bin/env nextflow

/*
 * created by Johanna Wong 6 August 2020
 */

/*
Parameter
/*
//run table files with file name and the file to read1 and read2*/
params.index = "$params.InputDir/index.csv"

//path to reference genome
params.shcontig = "$params.InputDir/reference/SH_assembly.fasta"
params.wgcontig = "$params.InputDir/reference/WG_assembly.fasta"
shcontig = file(params.shcontig)
wgcontig = file(params.wgcontig)


//setting for read trimming
params.mean_quality = 15
params.trimming_quality = 20

// Channels

Channel.fromPath(params.index)
        .splitCsv(header:true, sep: ',')
        .map{row-> tuple(row.sampleId, file(row.read1), file(row.read2))}
        .set { trim_ch }

//Step 1 - read trimming
  process fastp {
        tag "$sampleId"
        publishDir "${params.OutputDir}/trimmed", mode: 'copy'

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

// Step 2 - Read mapping to contigs
process mapping{

        tag "$sampleId"
        publishDir "${params.OutputDir}/bamfiles", mode: 'copy'

    input:
        set val(sampleId), file(tread1) from tread1_ch
        set val(sampleId), file(tread2) from tread2_ch

    output:
        file("*.bam")  into results2

    script:
        """
        bwa mem -t 16 $shcontig $tread1 $tread2 | samtools view -bS - > ${sampleId}_map2SH.bam
        bwa mem -t 16 $wgcontig $tread1 $tread2 | samtools view -bS - > ${sampleId}_map2WG.bam
        """
}

/* Step 2 - Read mapping to contigs
process mapping_WG{
        tag "$sampleId"
        conda 'bioconda::bwa bioconda::samtools openssl=1.0'
        publishDir "${params.OutputDir}/map2WG", mode: 'copy'

    input:
        set val(sampleId), file(tread1) from tread1_wg
        set val(sampleId), file(tread2) from tread2_wg

    output:
        set file("*.txt"), file("*.bam")  into results

    when: mapping for WG samples only
        sampleId =~ /EW[1-5]/ && sampleId != "EW10"

    script:
        """
        bwa mem $wgcontig $tread1 $tread2 | samtools view -S -b - > ${sampleId}_map2WG.bam
        samtools flagstat ${sampleId}_map2WG.bam > ${sampleId}_flagstat.txt
        """
}
*/
