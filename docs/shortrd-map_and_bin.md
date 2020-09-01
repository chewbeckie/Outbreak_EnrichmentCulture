# shortrd-map_and_bin.nf usage

Short reads are mapped to the assembled contigs (SH and WG) by bwa. The BAM files are then used for binning with metaBAT2.

1. Input the path to each individual illumina short readsets in `shortreads_*.csv`. 

2. Input the path to the assemblied genome using the `assembly-metaflye_and_polish.nf` script as `params.contig`. Also input a name for the name of the contig (`params.contigname`)

3. *(optional)* specify directory locations for input folder and output folder in `shortrd-map_and_bin.nf`

4. run the nextflow workflow script by using the following command (activate pbs if running on the UTS HPC)
```
nextflow run shortrd-map_and_bin.nf -c shortrd-map_and_bin.config -profile conda,pbs
```

## Dependence
* nextflow
* conda
    * bwa
    * fastp
    * samtools v1.9
    * metabat2