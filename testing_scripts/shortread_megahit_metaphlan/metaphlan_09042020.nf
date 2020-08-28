
Channel
  .fromFilePairs("/shared/homes/149306/metaphlan2_testing/input/*_{R1,R2}_001.fastq.gz")
  .set{trimmed_reads}

process metaphlan2 {

  conda '/shared/homes/149306/metaphlan2_testing/meta_env.yaml'
  cpus 8
  memory '12 GB'
  queue 'workq'
  executor 'pbspro'

      tag "$sampleId"
      publishDir "metaphlan2_out", mode: 'copy'

    input:
        set sampleId, file(reads) from trimmed_reads

    output:
        set val(sampleId), file("${sampleId}.bt2out") into metaphlan_res
        set val(sampleId), file("${sampleId}.mph2")

    script:

        """
        /shared/homes/149306/bin/metaphlan/metaphlan.py $reads --nproc 8 --input_type multifastq --bt2_ps very-sensitive \
        --bowtie2db /shared/homes/149306/bin/metaphlan/bowtie2db/mpa --bowtie2out ${sampleId}.bt2out > ${sampleId}.mph2
        """

}
