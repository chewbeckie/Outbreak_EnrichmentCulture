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
params.contig = "$params.InputDir/reference_genomes/SH_assembly_edited.fa" // "$params.InputDir/reference_genomes/WG_assembly_edited.fa"
params.contigname = "SH" //"WG"
contigs = file(params.contig)
location = params.contigname


//setting for read trimming
params.mean_quality = 15
params.trimming_quality = 20

// Channels
// Reads
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

// Step 2a - Read mapping to contigs
  process mapping {
        tag "$sampleId"
        publishDir "${params.OutputDir}/bamfiles", mode: 'copy'

    input:
        set val(sampleId), file(tread1) from tread1
        set val(sampleId), file(tread2) from tread2

    output:
        file("*map2${location}.bam")  into bam
        file("*.sorted.bam")  into (sortedbam, bam4vcf)

    script:
        """
        bwa mem -t 16 $contigs $tread1 $tread2 | \
         samtools view -bS - > ${sampleId}_map2${location}.bam
        samtools sort -@ 16 -o ${sampleId}_map2${location}.sorted.bam ${sampleId}_map2${location}.bam
        samtools index ${sampleId}_map2${location}.sorted.bam
        """
}

// Step 3 - metabat binning
process MetaBat2 {

        tag "$location"
        publishDir "${params.OutputDir}/metabat2", mode: 'copy'

    input:
        path '*'  from sortedbam.toList()

    output:
        file("metabat*")

    script:
        """
        jgi_summarize_bam_contig_depths --outputDepth Depth.txt --showDepth *
        metabat2 -i $contigs -a Depth.txt -o metabat_$location/bin
        mv Depth.txt metabat_$location
        """

}

// Step 4a - variant calling 
process varcall {
        conda "bioconda::lofreq" //lofreq has different requirment from other packages
        tag "$bam.simpleName"
        publishDir "${params.OutputDir}/varcall", mode: 'copy'

    input:
        path(bam) from bam4vcf

    output:
        file("*.vcf")

    script:
        """
        lofreq call -f $contigs -m 20 --no-default-filter -o ${bam.simpleName}.vcf $bam
        lofreq filter -v 3 -V 500 -i ${bam.simpleName}.vcf -o ${bam.simpleName}.filt.vcf
        """       
}

