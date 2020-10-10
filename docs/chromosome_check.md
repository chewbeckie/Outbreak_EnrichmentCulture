# chromosome_check.nf usage

The purpose of this script is to search the contigs against a database of bacterial chromosome genomes downloaded from NCBI.

User should supply the following:

* `genome_list` - a `txt` file containing the accession numbers of complete chromosome genomes from NCBI. This list would be used for creating the mash sketch for mash screen. Each entry is separated by a new line
* `query` - one or more `fasta` file(s) of contigs as query of mash search

The nextflow workflow script can be ran by using the following command (activate pbs if running on the UTS HPC)
```
chromosome_check.nf -c chromosome_check.config \ 
    --genome_list accession_no.txt --query contigs.fasta \
    -profile pbs
```

## Dependence
* nextflow
* conda
    * entrez-direct
    * mash