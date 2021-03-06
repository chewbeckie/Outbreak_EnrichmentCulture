#!/bin/bash

# Author: Johanna Wong

##################
# Set PBS Commands 
# PBS commands must come at the top of this script, before any other commands.
#################

# Set the resource requirements; 1 CPU, 5 GB memory and 5 minutes wall time.
#PBS -N bwa_SH
#PBS -l ncpus=28
#PBS -l mem=150GB
#PBS -l walltime=120:00:00 
#PBS -l centos6_node=yes

# There are several queues e.g. workq, smallq and others
#PBS -q workq

# Send email on abort, begin and end. 
# CHANGE 999777 to your staff or student number!
#PBS -m abe 
#PBS -M 149306@uts.edu.au

# define input
contigs="locA_assembly_edited.fa"
reads="interleaved_HiC_EW.fq"
out="HiC_EW_SHctg.bam"
work_dir="bin3c_Aug2020/"
# commands
source activate bwa
cd $work_dir
bwa mem -5SP -t 24 $contigs $reads | samtools view -F 0x904 -bS - > $out