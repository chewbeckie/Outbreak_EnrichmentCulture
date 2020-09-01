#!/bin/bash

# Author: Johanna Wong

##################
# Set PBS Commands 
# PBS commands must come at the top of this script, before any other commands.
#################

# Set the resource requirements; 1 CPU, 5 GB memory and 5 minutes wall time.
#PBS -N checkm
#PBS -l ncpus=8
#PBS -l mem=120GB
#PBS -l walltime=120:00:00 

# There are several queues e.g. workq, smallq and others
#PBS -q workq

# Send email on abort, begin and end. 
# CHANGE 999777 to your staff or student number!
#PBS -m abe 
#PBS -M 149306@uts.edu.au

work_dir="/shared/homes/149306/outbreak/enrichment_culture/work_in_progress/checkm_testing"

source activate checkm #version of checkm<1.1 
cd $work_dir
for i in input_files/*;
do prefix="$(basename -- $i)";
checkm lineage_wf -t 8 $i output_files/checkm_$prefix
checkm bin_qa_plot --image_type svg output_files/checkm_$prefix $i qa_$prefix
done
