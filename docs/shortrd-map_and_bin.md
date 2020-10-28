# shortrd-map_and_bin.nf usage

Short reads are mapped to the assembled contigs (locA and locB) by bwa. The BAM files are then used for binning with metaBAT2 and variant calling. This script requires the following as inputs:

* `--index` - a `.csv` runtable with information about the sample name, path to read1 and path to read 2 of each individual illumina short readsets
    ```
    sampleId,read1,read2
    sample1,path/to/r1,path/to/r2
    ```
* `--contig` - a path to the assembly fasta file
* `--contigname` - a user-defined name for the assembly
* `--out` - location for output files to be stored

The nextflow workflow script can be ran by using the following command (activate pbs if running on the UTS HPC)
```
shortrd-map_and_bin.nf -c shortrd-map_and_bin.config \ 
    --index path/to/runtable.csv --contig path/to/assembly_edited.fa --contigname locA \
    --out output_files -profile conda,pbs
```

## Dependence
* nextflow
* conda
    * bwa
    * fastp
    * samtools v1.9
    * metabat2
    * lofreq