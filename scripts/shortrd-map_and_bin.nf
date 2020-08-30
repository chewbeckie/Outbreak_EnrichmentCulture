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
params.shcontig = "$params.InputDir/reference_genomes/SH_assembly_edited.fa"
params.wgcontig = "$params.InputDir/reference_genomes/WG_assembly_edited.fa"
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
          set val(sampleId), file("trimmed_${sampleId}_R1.fastq.gz") into tread1
          set val(sampleId), file("trimmed_${sampleId}_R2.fastq.gz") into tread2
          file("fastp.*") into fastpresult

        script:
          """
          fastp -i "${read1}" -I "${read2}"\
            -o "trimmed_${sampleId}_R1.fastq.gz"\
            -O "trimmed_${sampleId}_R2.fastq.gz"
          """
}

// Splitting output into SH and WG groups
tread1.into{ tread1a; tread1b}
tread2.into{ tread2a; tread2b}

tread1a.filter{ it =~/trimmed_EW[1-5]/ && !(it =~ /trimmed_EW10/) }
      .set{tread1_wg}
tread2a.filter{ it =~/trimmed_EW[1-5]/ && !(it =~ /trimmed_EW10/) }
      .set{tread2_wg}

tread1b.filter{ (it =~/trimmed_EW[10, 6-9]/) && !(it =~ /trimmed_EW1_/)}
      .set{tread1_sh}
tread2b.filter{ (it =~/trimmed_EW[10, 6-9]/) && !(it =~ /trimmed_EW1_/)}
      .set{tread2_sh}


// Step 2a - Read mapping to contigs
process mapping_WG{
        tag "$sampleId"
        publishDir "${params.OutputDir}/bam_files", mode: 'copy'

    input:
        set val(sampleId), file(tread1) from tread1_wg
        set val(sampleId), file(tread2) from tread2_wg

    output:
        set val(sampleId), file("${sampleId}_map2WG.bam")  into bam_wg
        set val(sampleId), file("sorted_${sampleId}_map2WG.bam")  into sortedbam_wg

    //when: //mapping for WG samples only
    //    sampleId =~ /EW[1-5]/ && sampleId != "EW10"

    script:
        """
        bwa mem -t 16 $params.wgcontig $tread1 $tread2 | samtools view -bS - > ${sampleId}_map2WG.bam
        samtools sort -@ 16 -o sorted_${sampleId}_map2WG.bam ${sampleId}_map2WG.bam 
        """
}


// Step 2b - Read mapping to SH contigs
process mapping_SH{

        tag "$sampleId"
        publishDir "${params.OutputDir}/bam_files", mode: 'copy'

    input:
        set val(sampleId), file(tread1) from tread1_sh
        set val(sampleId), file(tread2) from tread2_sh

    output:
        set val(sampleId), file("${sampleId}_map2SH.bam")  into bam_sh
        set val(sampleId), file("sorted_${sampleId}_map2SH.bam")  into sortedbam_sh

    //when: //mapping for SH samples only
    //    sampleId =~ /EW[6-9]/ && sampleId == "EW10"

    script:
        """
        bwa mem -t 16 $params.shcontig $tread1 $tread2 | samtools view -bS - > ${sampleId}_map2SH.bam
        samtools sort -@ 16 -o sorted_${sampleId}_map2SH.bam ${sampleId}_map2SH.bam 
        """
}

// merge SH and WG channels back together for further processes
sortedbam_sh.mix(sortedbam_wg)
            .set{sortedbam}

sortedbam.view()
