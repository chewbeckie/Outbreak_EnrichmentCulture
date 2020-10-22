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

//Step 1 - read trimming
process fastp {
  conda 'bioconda::fastp'
  tag "$sampleId"
  publishDir "${params.out}/trimmed", mode: 'copy'
  memory '30 GB'
  cpus 12

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

//step 2: reformat input into interleaved fastq
process reformat {
    conda 'bioconda::bbmap'
    tag "$sampleId"
    publishDir "${params.out}/paired_reads", mode: 'copy'
    memory '40 GB'
    cpus 8

    input:
    set sampleId, file(read1), file(read2) from trimmed_reads

    output:
    set val(sampleId), file("interleaved_${sampleId}.fq") into paired_reads

    script:
    """
    reformat.sh in1=$read1 in2=$read2 out=interleaved_${sampleId}.fq
    """
}

//step 3: initial mapping for dedup
process init_map {
    conda 'bioconda::samtools bioconda::bwa'
    tag "$sampleId"
    publishDir "${params.out}/init_map", mode: 'copy'
    cpus 10
    memory '100 GB'
    input:
    set sampleId, file(reads) from paired_reads

    output:
    set val(sampleId), file("*.bam") into inimap_ch

    script:
    """
    # assume contigs are already indexed -- avoids different init_map processes
    # overwriting each other's index.
    bwa mem -t 9 $contigs -p $reads \
    | samtools view -bS - \
    | samtools sort -@ 1 -n - withdup_${sampleId}_map2${location}
    """
}


//step 4: dedup and bam2fq
// TODO: the cpu count could be abstracted to a variable here (and elsewhere)
process dedup {
    container 'docker://staphb/samtools:1.11'
    tag "$sampleId"
    publishDir "${params.out}/dedup_outputs", mode: 'copy'
    cpus 8
    memory '10 GB'
    input:
    set val(sampleId), file(withdup) from inimap_ch

    output:
    file("*.bam")
    set val(sampleId), file("*R1.fq"), file("*R2.fq") into dedup_ch

    """
    samtools fixmate -@ 7 -m $withdup - \
    | samtools sort -@ 7 - -T tmp0 -o tmp1.bam
    samtools markdup -@ 7 -r -s tmp1.bam tmp_dedup_${sampleId}_map2${location}.bam
    samtools sort -n -@ 8 tmp_dedup_${sampleId}_map2${location}.bam -o dedup_${sampleId}_map2${location}.namesort.bam
    samtools bam2fq dedup_${sampleId}_map2${location}.namesort.bam \
        -1 dedup_${sampleId}_map2${location}_R1.fq \
        -2 dedup_${sampleId}_map2${location}_R2.fq \
        -s singleton.fq
    rm -f tmp*.bam singleton.fq
    """
}


//step 5: map deduped reads to longread-assembled contigs

process map_reads {
    conda 'bioconda::samtools bioconda::bwa'

    tag "$sampleId"
    publishDir "${params.out}/bamfiles", mode: 'copy'
    cpus 9
    memory '140 GB'

    input:
    set val(sampleId), file(r1), file(r2) from dedup_ch

    output:
    file("*map2${location}.bam") into results
    file("*sorted.bam") into bam_ch

    script:
    """
    bwa mem -t 8 -5SP $contigs $r1 $r2 | samtools view -F 0x904 -bS - > ${sampleId}_map2${location}.bam
    samtools sort -@ 8 -n ${sampleId}_map2${location}.bam ${sampleId}_map2${location}.sorted
    rm -f tmp*
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
