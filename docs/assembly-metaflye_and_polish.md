# assembly-metaflye_and_polish.nf usage

De novo metagenome assembly using ONT long-read data and correction with illumina short-read data.

1. The individual illumina short-read libraries for each sample are merged by location using `gcat`
- WG (EW1-5)
- SH (EW6-10)

2. provide path to the ONT long read files in `longreads.csv` and the merged short read-sets in `shortreads.csv`, which are both located in the `input_files` folder.

3. setting the genome size (`params.gs`) in the `assembly-metaflye_and_polish.nf`, the genome size is set at 100m for the generation of WG/SH assemblies.

4. *(optional)* specify directory locations for input folder and output folder in `assembly-metaflye_and_polish.config`.

5. run the nextflow workflow script by using the following command (activate pbs if running on the UTS HPC).
```
nextflow run assembly-metaflye_and_polish.nf -c assembly-metaflye_and_polish.config -profile conda,pbs
```

## Dependence
* nextflow
* conda
* python