#!/usr/bin/env nextflow

/* created by Johanna Wong on 9 Jun 2020 */

//define parameters
//run table files with file name and the path to read1 and read2*/
params.index = "$baseDir/runtable.csv"
params.outdir = "$baseDir/output_files"
params.shcontigs = "locA_assembly_edited.fa"
params.wgcontigs = "locB_assembly_edited.fa"

//setting for read trimming
params.mean_quality = 15
params.trimming_quality = 20

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
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'


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
  publishDir "${params.outdir}/trimmed", mode: 'copy'

  input:
    set sampleId, path(read1), path(read2) from trim_ch
    val qual from params.mean_quality
    val trim_qual from params.trimming_quality

  output:
    set val(sampleId), path("trimmed_${sampleId}_R1.fastq.gz"), path("trimmed_${sampleId}_R2.fastq.gz") into trimmed_reads
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

    conda 'agbiome::bbtools'
    tag "$sampleId"
    publishDir "${params.outdir}/paired_reads", mode: 'copy'
    
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
    queue 'workq'
    memory '100 GB'
    cpus 20
    executor 'pbs'

    conda '/shared/homes/149306/miniconda3/envs/bwa'
    tag "$sampleId"
    publishDir "${params.outdir}/mapped_bams", mode: 'copy'


    input:
    set sampleId, path(reads) from paired_reads

    output:
    file("*.bam") into results

    script:
    """
    bwa mem -5SP $params.wgcontigs $reads | \
    samtools view -F 0x904 -bS - | \
    samtools sort -n -o ${sampleId}_WGctg.bam
    bwa mem -5SP $params.shcontigs $reads | \
    samtools view -F 0x904 -bS - | \
    samtools sort -n -o ${sampleId}_SHctg.bam
    """

}