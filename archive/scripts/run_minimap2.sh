#!/bin/bash
contig="SH_assembly.fasta"
reads="HS.pass.fastq"

#minimap2 -ax map-ont $contig $reads > promethion2SH.sam
samtools view -S -b promethion2SH.sam > promethion2SH.bam
samtools sort promethion2SH.bam -o promethion2SH.sorted.bam
samtools index promethion2SH.sorted.bam