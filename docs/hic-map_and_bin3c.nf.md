# hic-map_and_bin3c.nf usage

HiC readset are first preprocessed and merged into interleaved format and then mapped to assembled contigs (SH and WG) by bwa. The bam files are then used for variant calling.

1. Prepare a `.csv` runtable with information about the sample name, path to read1 and path to read 2 of each individual HiC readsets

```
sampleId,read1,read2
hiC_1,path/to/r1,path/to/r2
```

2. *(optional)* specify directory locations for input folder and output folder in `hic-map_and_bin3c.nf`

3. run the nextflow workflow script by using the following command (activate pbs if running on the UTS HPC)
```
hic-map_and_bin3c.nf -c hic-map_and_bin3c.config \
 --index input_files/hic_index.csv --contig input_files/SH_assembly_edited.fa --contigname SH -profile conda,pbs
```

## Dependence
* nextflow
* conda
    * bwa
    * fastp
    * samtools v1.9