#!/bin/bash

# Author: Johanna Wong

##################
# Set PBS Commands 
# PBS commands must come at the top of this script, before any other commands.
#################

# Set the resource requirements; 1 CPU, 5 GB memory and 5 minutes wall time.
#PBS -N metabat2
#PBS -l ncpus=16
#PBS -l mem=50GB
#PBS -l walltime=120:00:00 
#PBS -l centos6_node=yes

# There are several queues e.g. workq, smallq and others
#PBS -q medq

# Send email on abort, begin and end. 
# CHANGE 999777 to your staff or student number!
#PBS -m abe 
#PBS -M 149306@uts.edu.au

#set variables
workdir="/shared/homes/s1/outbreak/enrichment_culture/work_in_progress/metabat2_testing"
contigs="/shared/homes/s1/outbreak/enrichment_culture/aug2020-assembilies/SH_assembly_edited.fa"

#command
cd $workdir
source activate metabat2
jgi_summarize_bam_contig_depths --outputDepth Depth.txt --showDepth $workdir/sorted_bamfiles/sorted_*_map2SH.bam
metabat2 -i $contigs -a Depth.txt -o bins_dir/bin