# contigs_annotate.nf

The purpose of this script is to characterize the contigs assembled from `assembly_map_vcall.nf`. Taxonomic classification is performed with kraken. And the contigs were searched against PlsDB and a list of complete Enterobacteriacae genomes downloaded from ncbi (`data/input_files/complete_assembly_enterobacteriaceae_chromosome_info.csv`). Finally, AMR and pMLST gene screening are performed using ABRicate.

## Important to-do list before using the script

The following would need to be done before running the script:
* download and uncompress the kraken2 database (download [here](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads). Recommand using the Minikraken2 database as it is smaller)
* download fasta files of the chromosome/plasmid genomes listed on `complete_assembly_enterobacteriaceae_chromosome_info.csv` and `plsdb_20200629.csv`

## Usage
User should supply the following (please don't use relative path, it won't work):

* `name` - user-defined contigs name for output
* `contigs` - path to fasta file of the contigs
* `krakendb` - path to kraken database
* `plasmidlist` - csv file with path to the fasta files of the plasmid genome database collection
* `chromosomelist` - csv file with path to the fasta files of the chromosome genome database collection
* `contigs_as_ref` (optional) - if used, a 'reverse' mash screen would be performed using contigs as reference

The nextflow workflow script can be ran by using the following command 
```
./contigs_annotation.nf --name SH --contigs SH_assembly_edited.fa \
--krakendb minikraken2_v1_8GB --plasmidlist plasmids_list.tsv \ --chromosomelist input_files/chr_list.tsv 
```

## Dependence
* nextflow
* conda