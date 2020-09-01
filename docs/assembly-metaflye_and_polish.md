# assembly-metaflye_and_polish.nf usage

De novo metagenome assembly using ONT long-read data and correction with illumina short-read data.

1. provide path to the ONT long read files in `longreads_*.csv` and the merged short read-sets in `shortreads_*.csv`, which are both located in the `input_files` folder.

2. set the genome size (`params.gs`) and contig name (`params.contigname`) in the `assembly-metaflye_and_polish.nf`, the genome size is set at 100m for the generation of WG/SH assemblies.

4. *(optional)* specify directory locations for input folder and output folder in `assembly-metaflye_and_polish.config`.

5. run the nextflow workflow script by using the following command (activate pbs if running on the UTS HPC).
```
nextflow run assembly-metaflye_and_polish.nf -c assembly-metaflye_and_polish.config -profile conda,pbs
```

## Dependence
* nextflow
* conda (dependencies are listed in `assembly-metaflye_and_polish.config`)
    * flye=2.7
    * fastp
    * ntedit=1.3.2
* python