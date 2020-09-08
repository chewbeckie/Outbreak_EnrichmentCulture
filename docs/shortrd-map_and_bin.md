# shortrd-map_and_bin.nf usage

Short reads are mapped to the assembled contigs (SH and WG) by bwa. The BAM files are then used for binning with metaBAT2 and variant calling.

1. Prepare a `.csv` runtable with information about the sample name, path to read1 and path to read 2 of each individual illumina short readsets

```
sampleId,read1,read2
sample1,path/to/r1,path/to/r2
```

2. *(optional)* specify directory locations for input folder and output folder in `shortrd-map_and_bin.nf`

3. run the nextflow workflow script by using the following command (activate pbs if running on the UTS HPC)
```
shortrd-map_and_bin.nf -c shortrd-map_and_bin.config \ 
    --path/to/runtable.csv --contig path/to/assembly_edited.fa --contigname SH \
    -profile conda,pbs
```

## Dependence
* nextflow
* conda
    * bwa
    * fastp
    * samtools v1.9
    * metabat2