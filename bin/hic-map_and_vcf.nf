#!/usr/bin/env nextflow

/* created by Johanna Wong on 9 Jun 2020 */

/* define parameters*/
//output directory location
params.out
//run table files with file name and the path to read1 and read2 of HiC data*/
params.index //="$params.InputDir/hicreads.csv"
//path to reference genome
params.contig  //="$params.InputDir/reference_genomes/SH_assembly_edited.fa" or "$params.InputDir/reference_genomes/WG_assembly_edited.fa"
params.contigname //="SH" //"WG"
contigs = file(params.contig)
location = params.contigname

//Channels
Channel.fromPath(params.index)
        .splitCsv(header:true, sep: ',')
        .map{row-> tuple(row.sampleId, file(row.read1), file(row.read2))}
        .into { samples_ch; trim_ch }

//step 1: clean reads

/*process clean_reads{
    queue 'smallq'
    memory '24 GB'
    cpus 4
    executor 'pbs'

    conda 'agbiome::bbtools'
    tag "$sampleId"
    publishDir "${params.out}/trimmed_reads", mode: 'copy'


    input:
    set sampleId, file(read1), file(read2) from samples_ch

    output:
    set val(sampleId), file("clean_${sampleId}_1.fq"), file("clean_${sampleId}_2.fq") into trimmed_reads

    script:
    """
    bbduk.sh in1=$read1 in2=$read2 ref=$params.adapterseq\
    k=23 hdist=1 mink=11 ktrim=r tpe tbo \
    ftm=5 qtrim=r trimq=10 \
    threads=4 -Xmx20g\
    out1=clean_${sampleId}_1.fq out2=clean_${sampleId}_2.fq
    """
}
*/

//Step 2 - read trimming
process fastp {
  tag "$sampleId"
  publishDir "${params.out}/trimmed", mode: 'copy'

  input:
    set sampleId, path(read1), path(read2) from trim_ch

  output:
    set val(sampleId), path("*_R1.fastq.gz"), path("*_R2.fastq.gz") into trimmed_reads
    file("fastp.*")

  script:
    """
    fastp -i "${read1}" -I "${read2}"\
      -o "trimmed_${sampleId}_R1.fastq.gz"\
      -O "trimmed_${sampleId}_R2.fastq.gz"
    """
}

//step 3: reformat input into interleaved fastq

process reformat {

    tag "$sampleId"
    publishDir "${params.out}/paired_reads", mode: 'copy'
    
    input:
    set sampleId, file(read1), file(read2) from trimmed_reads

    output:
    set val(sampleId), file("interleaved_${sampleId}.fq") into paired_reads

    script:
    """
    reformat.sh in1=$read1 in2=$read2 out=interleaved_${sampleId}.fq
    """
}

//step 4: map reads to longread-assembled contigs

process map_reads {

    tag "$sampleId"
    publishDir "${params.out}/bamfiles", mode: 'copy'


    input:
    set sampleId, file(reads) from paired_reads

    output:
    file("*map2${location}.bam") into results
    file("*sorted.bam") into bam_ch

    script:
    """
    bwa mem -t 16 -5SP $contigs $reads | samtools view -F 0x904 -bS - > ${sampleId}_map2${location}.bam
    samtools sort -@ 16 -o ${sampleId}_map2${location}.sorted.bam -n ${sampleId}_map2${location}.bam
    """
}

// Step 5 - variant calling 
process varcall {
        conda "bioconda::lofreq" //lofreq has different requirment from other packages
        tag "$bam.simpleName"
        publishDir "${params.out}/varcall", mode: 'copy'

    input:
        path(bam) from bam_ch

    output:
        file("*.vcf") into vcf_ch

    script:
        """
        lofreq call -f $contigs -m 20 --no-default-filter -o ${bam.simpleName}.vcf $bam
        lofreq filter -v 3 -V 500 -i ${bam.simpleName}.vcf -o ${bam.simpleName}.filt.vcf
        """
}

// Step 5 - index and compress vcf
process vcf_index {
        tag "$vcf.simpleName"
        publishDir "${params.out}/varcall", mode: 'copy'
    
    input:
        file(vcf) from vcf_ch

    output:
        file("*")
    
    script:
        """
        bgzip --index -@ 16 $vcf
        """
}