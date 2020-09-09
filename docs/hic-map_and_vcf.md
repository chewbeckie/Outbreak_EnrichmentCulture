# hic-map_and_vcf.nf usage

HiC readset are first preprocessed and merged into interleaved format and then mapped to assembled contigs (SH and WG) by bwa. The bam files are then used for variant calling. This script requires the following as inputs:

*  `--index` - a `.csv` runtable with information about the sample name, path to read1 and path to read 2 of each individual HiC readsets
    ```
    sampleId,read1,read2
    hiC_1,path/to/r1,path/to/r2
    ```

* `--contig` - a path to the assembly fasta file
* `--contigname` - a user-defined name for the assembly
* `--out` - location for output files to be stored


The nextflow workflow script can be ran by using the following command (activate pbs if running on the UTS HPC).
```
hic-map_and_vcf.nf -c hic-map_and_vcf.config \
 --index input_files/hic_index.csv --contig input_files/SH_assembly_edited.fa --contigname SH \
 --out output_files -profile conda,pbs
```

## Dependence
* nextflow
* conda
    * bwa
    * fastp
    * samtools v1.9
    * lofreq
    * bbtools