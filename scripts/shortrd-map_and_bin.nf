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

//setting for read trimming
params.mean_quality = 15
params.trimming_quality = 20

// Channels
// Reads
Channel.fromPath(params.index)
        .splitCsv(header:true, sep: ',')
        .map{row-> tuple(row.sampleId, file(row.read1), file(row.read2))}
        .set { trim_ch }

//Contigs
Channel.from("SH", file(params.shcontig),
             "WG", file(params.wgcontig))
       .buffer( size: 2)
       .set {contig_ch}

//Step 1 - read trimming
  process fastp {
        tag "$sampleId"
        publishDir "${params.OutputDir}/trimmed", mode: 'copy'

        input:
          set sampleId, file(read1), file(read2) from trim_ch

        output:
          set val(sampleId), file("trimmed_${sampleId}_R1.fastq.gz") into tread1a, tread1b
          set val(sampleId), file("trimmed_${sampleId}_R2.fastq.gz") into tread2a, tread2b
          file("fastp.*") into fastpresult

        script:
          """
          fastp -i "${read1}" -I "${read2}"\
            -o "trimmed_${sampleId}_R1.fastq.gz"\
            -O "trimmed_${sampleId}_R2.fastq.gz"
          """
}

// Splitting output into SH and WG groups
tread1a.filter{ it =~/trimmed_EW[1-5]/ && !(it =~ /trimmed_EW10/) }
      .set{tread1_wg}
tread2a.filter{ it =~/trimmed_EW[1-5]/ && !(it =~ /trimmed_EW10/) }
      .set{tread2_wg}

tread1b.filter{ (it =~/trimmed_EW[10, 6-9]/) && !(it =~ /trimmed_EW1_/)}
      .set{tread1_sh}
tread2b.filter{ (it =~/trimmed_EW[10, 6-9]/) && !(it =~ /trimmed_EW1_/)}
      .set{tread2_sh}


// Step 2a - Read mapping to WG contigs
  process mapping_WG {
        tag "$sampleId"
        publishDir "${params.OutputDir}/bamfiles", mode: 'copy'

    input:
        set val(sampleId), file(tread1) from tread1_wg
        set val(sampleId), file(tread2) from tread2_wg

    output:
        file("${sampleId}_map2WG.bam")  into bam_wg
        file("sorted_${sampleId}_map2WG.bam")  into (sortedbam_wg, bam4vcf_wg)

    script:
        """
        bwa mem -t 16 $params.wgcontig $tread1 $tread2 | samtools view -bS - > ${sampleId}_map2WG.bam
        samtools sort -@ 16 -o sorted_${sampleId}_map2WG.bam ${sampleId}_map2WG.bam
        samtools index sorted_${sampleId}_map2WG.bam
        """
}


// Step 2b - Read mapping to SH contigs
process mapping_SH{

        tag "$sampleId"
        publishDir "${params.OutputDir}/bamfiles", mode: 'copy'

    input:
        set val(sampleId), file(tread1) from tread1_sh
        set val(sampleId), file(tread2) from tread2_sh

    output:
        file("${sampleId}_map2SH.bam")  into bam_sh
        file("sorted_${sampleId}_map2SH.bam")  into (sortedbam_sh, bam4vcf_sh)

    script:
        """
        bwa mem -t 16 $params.shcontig $tread1 $tread2 | samtools view -bS - > ${sampleId}_map2SH.bam
        samtools sort -@ 16 -o sorted_${sampleId}_map2SH.bam ${sampleId}_map2SH.bam
        samtools index sorted_${sampleId}_map2SH.bam
        """
}


// concatenate bam output from both SH and WG set
sortedbam_sh.concat(sortedbam_wg)
            .into{ sortedbam_all ; sortedbam }

// Step 3 - metabat binning
process MetaBat2 {

        tag "$location"
        publishDir "${params.OutputDir}/metabat2", mode: 'copy'

    input:
        path '*'  from sortedbam_all.toList()
        set val(location), file(contig) from contig_ch

    output:
        file("metabat*")

    script:
        """
        jgi_summarize_bam_contig_depths --outputDepth Depth.txt --showDepth *_map2${location}.bam
        metabat2 -i $contig -a Depth.txt -o metabat_$location/bin
        mv Depth.txt metabat_$location
        """

}

// Step 4a - variant calling for WG samples
process varcall_WG {
        conda "bioconda::lofreq" //lofreq has different requirment from other packages
        tag "$bam.simpleName"
        publishDir "${params.OutputDir}/varcall", mode: 'copy'

    input:
        path(bam) from bam4vcf_wg

    output:
        file("*.vcf")
        file("sorted*")

    script:
        """
        lofreq call -f $params.wgcontig -m 20 --no-default-filter -o ${bam.simpleName}.vcf $bam
        lofreq filter -v 3 -V 500 -i ${bam.simpleName}.vcf -o ${bam.simpleName}.filt.vcf
        """       
}

// Step 4b - variant calling for SH samples
process varcall_SH {
        conda "bioconda::lofreq" //lofreq has different requirment from other packages
        tag "$bam.simpleName"
        publishDir "${params.OutputDir}/varcall", mode: 'copy'

    input:
        path(bam) from bam4vcf_sh

    output:
        file("*.vcf")
        file("sorted*")

    script:
        """
        lofreq call -f $params.shcontig -m 20 --no-default-filter -o ${bam.simpleName}.vcf $bam
        lofreq filter -v 3 -V 500 -i ${bam.simpleName}.vcf -o ${bam.simpleName}.filt.vcf
        """       
}
