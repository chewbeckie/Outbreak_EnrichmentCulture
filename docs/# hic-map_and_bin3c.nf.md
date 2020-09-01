# hic-map_and_bin3c.nf usage

HiC readset are first preprocessed and merged into interleaved format and then mapped to assembled contigs (SH and WG) by bwa.

1. Input the path to each individual HiC readsets in `hicreads.csv`. 

2. Input the path to the assemblied genome using the `hic-map_and_bin3c.nf` script as `params.contig`. Also input a name for the name of the contig (`params.contigname`)

3. *(optional)* specify directory locations for input folder and output folder in `hic-map_and_bin3c.nf`

4. run the nextflow workflow script by using the following command (activate pbs if running on the UTS HPC)
```
nextflow run hic-map_and_bin3c.nf -c hic-map_and_bin3c.config -profile conda,pbs
```

## Dependence
* nextflow
* conda
    * bwa
    * fastp
    * samtools v1.9