#!/bin/bash
#PBS -l ncpus=12
#PBS -l walltime=36:00:00
#PBS -l mem=60g
#PBS -N varcall

# There are several queues e.g. workq, smallq and others
#PBS -q medq

# Send email on abort, begin and end. 
# CHANGE 999777 to your staff or student number!
#PBS -m abe 
#PBS -M 149306@uts.edu.au

work_dir="/shared/homes/149306/outbreak/enrichment_culture/work_in_progress/varcall_testing"
ref="/shared/homes/s1/outbreak/enrichment_culture/aug2020-assembilies/SH_assembly_edited.fa"

source activate lofreq
cd $work_dir
for bam in $work_dir/bamfiles_shortreads/sorted*_map2SH.bam;
do bname=`basename $bam .bam`
samtools index $bam
lofreq call -f $ref -m 20 --no-default-filter -o $bname.vcf $bam
lofreq filter -v 3 -V 500 -i $bname.vcf -o $bname.filt.vcf
done