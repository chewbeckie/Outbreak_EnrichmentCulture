# contigs_annotate.nf

The purpose of this script is to characterize the contigs assembled from `assembly_map_vcall.nf`. Taxonomic classification is performed with kraken. And the contigs were searched against PlsDB and a list of complete Enterobacteriacae genomes downloaded from ncbi (`data/input_files/complete_assembly_enterobacteriaceae_chromosome_info.csv`). Finally, AMR and pMLST gene screening are performed using ABRicate.

## Important to-do list before using the script

The following would need to be done before running the script:
* download and uncompress the kraken2 database (download [here](https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads). Recommand using the Minikraken2 database as it is smaller)
* download fasta files of the chromosome/plasmid genomes listed on `complete_assembly_enterobacteriaceae_chromosome_info.csv` and `plsdb.tsv`. These files can be downloaded by using the `dl_fasta.nf` script in bin

## Usage
User should supply the following (please don't use relative path, it won't work):

* `name` - user-defined contigs name for output
* `contigs` - path to fasta file of the contigs
* `krakendb` - path to kraken database
* `plspath` - path + matching glob pattern of the plasmid genomes database fasta files in the directory
* `chrpath` - path + matching glob pattern of the chromosome genomes database fasta files in the directory
* `contigs_as_ref` (optional) - if used, a 'reverse' mash screen would be performed using contigs as reference instead of query
* `chr_info` - path to chromosome information file e.g. complete_assembly_enterobacteriaceae_chromosome_info.csv
* `pls_info` - path to plasmid information file e.g. plsdb.tsv

* `plsmsh` - path to plasmid database .msh mash sketch file, if supplied, `contigs_as_ref` and `plspath` must be absent
* `chrmsh` - path to chromosome database .msh mash sketch file, if supplied, `contigs_as_ref` and `chrpath` must be absent

The nextflow workflow script can be ran by using the following command 
```
./contigs_annotation.nf --name locA --contigs locA_assembly_edited.fa --krakendb minikraken2_v1_8GB \
--plspath input_files/plsdb_genomes/*.fa --chrpath input_files/chr_genomes/*.fa \
--chr_info complete_assembly_enterobacteriaceae_chromosome_info.csv --pls_ino plsdb.tsv
```

Alternatively, if .msh sketches are available, use this command
```
./contigs_annotation.nf --name locA --contigs locA_assembly_edited.fa --krakendb minikraken2_v1_8GB \
--plsmsh plsdb.msh --chrpath chr.msh \
--chr_info complete_assembly_enterobacteriaceae_chromosome_info.csv --pls_ino plsdb.tsv
```

## Dependence
* nextflow
* conda