#!/bin/bash

# Author: Johanna Wong

##################
# Set PBS Commands 
# PBS commands must come at the top of this script, before any other commands.
#################

# Set the resource requirements; 1 CPU, 5 GB memory and 5 minutes wall time.
#PBS -N 20200708_ntEdit_polish
#PBS -l ncpus=16
#PBS -l mem=200gb
#PBS -l walltime=12:00:00 
#PBS -l centos6_node=yes

# There are several queues e.g. workq, smallq and others
#PBS -q workq

# Send email on abort, begin and end. 
# CHANGE 999777 to your staff or student number!
#PBS -m abe 
#PBS -M 149306@uts.edu.au

###############
# setting variables
###############
r1="/shared/homes/149306/metaflye_testing/ntedit_polish/merged_SH_R1.fq.gz"
r2="/shared/homes/149306/metaflye_testing/ntedit_polish/merged_SH_R2.fq.gz"
contigs="/shared/homes/149306/metaflye_testing/SH_metaflye/assembly.fasta"
#bloomfilter=/shared/homes/149306/metaflye_testing/ntedit_polish/solidBF_k40_SH.bf

###############
# Start the Job
###############
source activate ntedit
# step 1 run nthits
nthits -c 1 --outbloom -p solidBF_k40_SH -b 36 -k 40 -t 16 $r1 $r2
# step 2 run ntedit
#ntedit -m 1 -f $contigs -r $bloomfilter -b ntEdit_SH_assembly