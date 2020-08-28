#!/usr/bin/env nextflow

/*
 * created by Johanna Wong 17 April 2020
 * wastewater samples - read mapping after coassembly
 */
params.outdir_sh = "/shared/homes/149306/megahit_testing/output_coassembly_SH"
params.outdir_wg = "/shared/homes/149306/megahit_testing/output_coassembly_WG"
contigs_sh = "/shared/homes/149306/megahit_testing/output_coassembly_SH/SH.contigs.fa"
contigs_wg = "/shared/homes/149306/megahit_testing/output_coassembly_WG/WG.contigs.fa"

// coassembly of EW 1-3 and EW 6-8
Channel
  .from(6..8)
  .map{ chr -> ["EW$chr",
                file("/shared/homes/149306/megahit_testing/input_coassembly/trimmed_EW${chr}_*_R1_001.fastq"),
                file("/shared/homes/149306/megahit_testing/input_coassembly/trimmed_EW${chr}_*_R2_001.fastq")]}
  .set{trimmed_sh}

Channel
  .from(1..3)
  .map{ chr -> ["EW$chr",
                file("/shared/homes/149306/megahit_testing/input_coassembly/trimmed_EW${chr}_*_R1_001.fastq"),
                file("/shared/homes/149306/megahit_testing/input_coassembly/trimmed_EW${chr}_*_R2_001.fastq")]}
  .set{trimmed_wg}

process mapping_sh{
  queue 'workq'
  memory '20 GB'
  cpus  8
  executor 'pbspro'

        tag "$sampleId"
        conda 'bioconda::bwa; bioconda::samtools'
        publishDir "${params.outdir_sh}/bwa", mode: 'copy'

    input:
        set val(sampleId), file(read1), file(read2) from trimmed_sh

    output:
        set file("*.txt"), file("*.bam")  into results_sh

    script:
        """
        bwa index $contigs_sh -p SH.contigs.fa
        bwa mem SH.contigs.fa $read1 $read2 | samtools view -S -b - > ${sampleId}_mapped.bam
        samtools flagstat ${sampleId}_mapped.bam > ${sampleId}_flagstat.txt
        """
}

process mapping_wg{
  queue 'workq'
  memory '20 GB'
  cpus  8
  executor 'pbspro'

        tag "$sampleId"
        conda 'bioconda::bwa; bioconda::samtools'
        publishDir "${params.outdir_wg}/bwa", mode: 'copy'

    input:
        set val(sampleId), file(read1), file(read2) from trimmed_wg

    output:
        set file("*.txt"), file("*.bam")  into results_wg

    script:
        """
        bwa index $contigs_wg -p WG.contigs.fa
        bwa mem WG.contigs.fa $read1 $read2 | samtools view -S -b - > ${sampleId}_mapped.bam
        samtools flagstat ${sampleId}_mapped.bam > ${sampleId}_flagstat.txt
        """
}
