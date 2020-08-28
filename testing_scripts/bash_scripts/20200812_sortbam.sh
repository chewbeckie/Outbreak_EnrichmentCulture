#!/bin/bash

# Author: Johanna Wong

##################
# Set PBS Commands 
# PBS commands must come at the top of this script, before any other commands.
#################

# Set the resource requirements; 1 CPU, 5 GB memory and 5 minutes wall time.
#PBS -N samsort
#PBS -l ncpus=16
#PBS -l mem=50GB
#PBS -l walltime=120:00:00 
#PBS -l centos6_node=no

# There are several queues e.g. workq, smallq and others
#PBS -q i3

# Send email on abort, begin and end. 
# CHANGE 999777 to your staff or student number!
#PBS -m abe 
#PBS -M 149306@uts.edu.au

#set variables
workdir="/shared/homes/149306/outbreak/enrichment_culture/work_in_progress/bwa_testing/output_files/bamfiles"

cd $workdir
source activate bwa
for i in $workdir/EW*.bam; 
do prefix=$(basename $i .bam);
samtools
samtools sort -@ 18 -o sorted_$prefix.bam $i ; done