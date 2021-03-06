# hic-map_and_vcf.nf usage

HiC readset are first preprocessed and merged into interleaved format and then mapped to assembled contigs (locA and locB) by bwa. Deduplication was then performed on the bam files, the deduped bam files were then converted back to fastq files. Then, the fastq files were remapped onto the assembled contigs and then used for variant calling. This script requires the following as inputs:

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
 --index input_files/hic_index.csv --contig input_files/locA_assembly_edited.fa --contigname locA \
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