#!/bin/bash

# Author: Johanna Wong

##################
# Set PBS Commands 
# PBS commands must come at the top of this script, before any other commands.
#################

# Set the resource requirements; 1 CPU, 5 GB memory and 5 minutes wall time.
#PBS -N 20200721_checkV
#PBS -l ncpus=20
#PBS -l mem=60gb
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
contigs="/shared/homes/s1/outbreak/enrichment_culture/metaflye/ntEdit_polish/ntEdit_SH_assembly_edited.fa"
output_directory="/shared/homes/149306/checkv_testing/SH_checkV_output"
database="/shared/homes/149306/checkv-db-v0.6"

###############
# Start the Job
###############
source activate checkv
checkv end_to_end $contigs $output_directory -t 16 -d $database