#!/bin/bash
# date = 20200511

# Author: Johanna Wong

##################
# Set PBS Commands 
# PBS commands must come at the top of this script, before any other commands.
#################

# Set the resource requirements; 1 CPU, 5 GB memory and 5 minutes wall time.
#PBS -N kraken_job
#PBS -l ncpus=5
#PBS -l mem=100GB
#PBS -l walltime=12:00:00 

# There are several queues e.g. workq, smallq and others
#PBS -q workq

# Send email on abort, begin and end. 
# CHANGE 999777 to your staff or student number!
#PBS -m abe 
#PBS -M 149306@uts.edu.au

# define variable
workdir="/shared/homes/149306/outbreak/enrichment_culture/work_in_progress/kraken_testing/SH_kraken2_out"
db="/shared/homes/149306/bin/kraken2_db"
contig="/shared/homes/149306/outbreak/enrichment_culture/work_in_progress/kraken_testing/SH_assembly_edited.fa"

# command
cd $workdir
kraken2 --db $db $contig --threads 4 \
--unclassified-out unclassified.tsv --classified-out classified.tsv \
--memory-mapping --report report.tsv --output output.tsv --use-names
