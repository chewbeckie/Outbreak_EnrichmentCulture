#!/bin/bash

# Author: Johanna Wong

##################
# Set PBS Commands 
# PBS commands must come at the top of this script, before any other commands.
#################

# Set the resource requirements; 1 CPU, 5 GB memory and 5 minutes wall time.
#PBS -N 20200617_metaflye_SH_job
#PBS -l ncpus=16
#PBS -l mem=40GB
#PBS -l walltime=200:00:00 

# There are several queues e.g. workq, smallq and others
#PBS -q medq

# Send email on abort, begin and end. 
# CHANGE 999777 to your staff or student number!
#PBS -m abe 
#PBS -M 149306@uts.edu.au

# define input
reads="/shared/homes/149306/metaflye_testing/promethion_run/HS.pass.fastq"
outdir="/shared/homes/149306/metaflye_testing/SH_metaflye"

# commands
source activate metaflye
flye --nano-raw $reads \
--meta --genome-size 5m \
-t 16 -o $outdir