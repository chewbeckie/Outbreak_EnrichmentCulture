#!/usr/bin/env nextflow

/*
 * created by Johanna Wong Sep 2020
 */

/*Parameters-----------------------------------------------------------------------------------------------------------------------------*/
// set output directory path
params.out //= "output_files"

// run table files with sample name and paths to read1, read2 
params.shortreads //= "$params.InputDir/shortreads_SH.csv" or "$params.InputDir/shortreads_WG.csv"
// sample name and path to long read files 
params.longreads //= "$params.InputDir/longreads_SH.csv" or "$params.InputDir/longreads_WG.csv"

// user-defined genome names
params.contigname //= "SH" or "WG"
// path to genome file (if supplied, assembly steps will be skipped)
params.contig //= "$params.InputDir/reference_genomes/SH_assembly_edited.fa" or "$params.InputDir/reference_genomes/WG_assembly_edited.fa"
// setting for metaflye
params.gs = "100m"

/*Condition validations*/
//if contigs were not provided, conditional assembly steps will run and the newly assembled contig will be used for contig for mapping
if(params.contig){ contig_flag = true } else { contig_flag = false }
if(contig_flag == true){contigs = file(params.contig)} else { contigs = file(newcontig) }

/*Channels-----------------------------------------------------------------------------------------------------------------------------*/

Channel.fromPath(params.shortreads)
        .splitCsv(header:true, sep: ',')
        .map{row-> tuple(row.sampleId, file(row.read1), file(row.read2))}
        .set { trim_ch }

Channel.fromPath(params.longreads)
        .splitCsv(header:true, sep: ',')
        .map{row-> tuple(row.sampleId, file(row.longread))}
        .set { longread_ch }


/*Processes-----------------------------------------------------------------------------------------------------------------------------*/

//Part 1: preprocessing

//Step 1.1 Read trimming
process fastp {
  tag "$sampleId"
  publishDir "${params.out}/trimmed_reads", mode: 'copy'

  input:
    tuple sampleId, file(read1), file(read2) from trim_ch

  output:

    file("trimmed_${sampleId}_R1.fastq.gz") into tread1_ch
    file("trimmed_${sampleId}_R2.fastq.gz") into tread2_ch
    tuple val(sampleId), file("trimmed_${sampleId}_R1.fastq.gz") into tr1_ch
    tuple val(sampleId), file("trimmed_${sampleId}_R2.fastq.gz") into tr2_ch

    file("*_fastp.*")

  script:
    """
    fastp -i "${read1}" -I "${read2}"\
      -o "trimmed_${sampleId}_R1.fastq.gz"\
      -O "trimmed_${sampleId}_R2.fastq.gz"\
      -j ${sampleId}_fastp.json -h ${sampleId}_fastp.html
    """
}

//Part 2 (Optional): genome assembly with longreads and polishing

//Step 2.1 Merging read1 and read2 together (only when contig not provided)
process merge_reads{
  publishDir "${params.out}/trimmed", mode: 'copy'

  input:
    path '*_R1.fastq.gz' from tread1_ch.toList()
    path '*_R2.fastq.gz' from tread2_ch.toList()
  
  output:
    tuple file("*_R1.fastq.gz"), file("*_R2.fastq.gz") into mergedreads_ch

  when:
    contig_flag == false

  script:
    """
    cat *_R1.fastq.gz > merged_${params.contigname}_R1.fastq.gz
    cat *_R2.fastq.gz > merged_${params.contigname}_R2.fastq.gz
    """

}

//Step 2.2 Metaflye assembly (only when contig not provided)
process metaflye_assembly{
  tag "$sampleId"
  publishDir "${params.out}/assembly", mode: 'copy'

  input:
    tuple val(sampleId), file(longread) from longread_ch

  output:
    tuple val(sampleId), file("flye_${sampleId}/*") into assembly_result
    tuple val(sampleId), file("${sampleId}_assembly.gfa") into graph_ch
  
  when:
    contig_flag == false

  script:
    """
    flye --nano-raw $longread \
      --meta --genome-size $params.gs \
      -t 16 -o flye_$sampleId
    cp flye_${sampleId}/assembly_graph.gfa ${sampleId}_assembly.gfa
    """
}

//Step 2.3 Polish with ntEdit (only when contig not provided)
process ntedit_polish{  
  tag "$sampleId"
  publishDir "${params.out}/polished_assembly", mode: 'copy'

  input:
    tuple file(merged_read1), file(merged_read2) from mergedreads_ch
    tuple val(sampleId), file(graph) from graph_ch

  output:
    file("*")
    file("${sampleId}_assembly_edited.fa") into newcontig

  when:
    contig_flag == false

  script:
    """
    flye_gfa_to_fasta_for_polishing.py $graph ${sampleId}_assembly.fa
    nthits -c 1 --outbloom -p solidBF_${sampleId} -b 36 -k 40 -t 16 $merged_read1 $merged_read2
    ntedit -m 1 -f ${sampleId}_assembly.fa -r solidBF_${sampleId}_k40.bf -b ${sampleId}_assembly
    flye_polished_fa_to_gfa.py $graph ${sampleId}_assembly_edited.fa ${sampleId}_assembly_edited.gfa
    """
}

// Part 3: Read mapping and variant calling

//Step 3.1 Read mapping to contigs
process mapping {
  tag "$sampleId"
  publishDir "${params.out}/bam_files", mode: 'copy'

  input:
    tuple val(sampleId), file(tread1) from tr1_ch
    tuple val(sampleId), file(tread2) from tr2_ch

  output:
    file("*map2${params.contigname}.bam")
    file("*.sorted.bam")  into (sortedbam, bam4vcf)

  script:
  """
    bwa mem -t 16 $contigs $tread1 $tread2 | \
    samtools view -bS - > ${sampleId}_map2${params.contigname}.bam
    samtools sort -@ 16 -o ${sampleId}_map2${params.contigname}.sorted.bam ${sampleId}_map2${params.contigname}.bam
    samtools index ${sampleId}_map2${params.contigname}.sorted.bam
  """
}

//Step 3.2 Metabat binning
process metaBat2 {
  tag "$params.contigname"
  publishDir "${params.out}/metabat2", mode: 'copy'

  input:
    file('*')  from sortedbam.toList()

  output:
    file("metabat*")

  script:
    """
    jgi_summarize_bam_contig_depths --outputDepth Depth.txt --showDepth *
    metabat2 -i $contigs -a Depth.txt -o metabat_$params.contigname/bin
    mv Depth.txt metabat_$params.contigname
    """
}

//Step 3.3 Variant calling 
process varcall {
  conda "bioconda::lofreq" //lofreq has different requirment from other packages
  tag "$bam.simpleName"
  publishDir "${params.out}/vcf_files", mode: 'copy'

  input:
    file(bam) from bam4vcf

  output:
    file("*.vcf") into vcf_ch

  script:
    """
    lofreq call -f $contigs -m 20 --no-default-filter -o ${bam.simpleName}.vcf $bam
    lofreq filter -v 3 -V 500 -i ${bam.simpleName}.vcf -o ${bam.simpleName}.filt.vcf
    """       
}

//Step 3.4 Index and compress vcf
process vcf_index {
  tag "$vcf.simpleName"
  publishDir "${params.out}/vcf_files", mode: 'copy'
    
  input:
    file(vcf) from vcf_ch

  output:
    file("*")
    
  script:
    """
    bgzip --index -@ 16 $vcf
    """
}