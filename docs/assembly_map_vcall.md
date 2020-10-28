# assembly_map_vcall.nf usage

This workflow is the integrated workflow for metagenomic genome assembly, read mapping, binning and variant calling. Firstly, steps from `assembly-metaflye_and_polish.nf` can be initiated optinally to prepare a genome assembly. Long reads were used to assembled into a draft genome using Flye, then short reads were used to polish the assembled genome. Afterwards, steps are same as `shortrd-map_and_bin.nf`. Short reads are mapped to the assembled contigs (locA and locB) by bwa. The BAM files are then used for binning with metaBAT2 and variant calling. This script requires the following as inputs:

* `--shortreads` - a `.csv` runtable with information about the sample name, path to read1 and path to read 2 of each individual illumina short readsets
    ```
    sampleId,read1,read2
    sample1,path/to/r1,path/to/r2
    ```
* `--longreads` - in same format as the shortreads runtable, a `.csv` runtable with information about the sample name, path to read1 and path to read 2 of each individual ONT long readsets
* `--contigname` - a user-defined name for the assembly
* `--contig` - (optional) a path to the assembly fasta file. When supplied, genome assembly and correction steps will be skipped.
* `--out` - location for output files to be stored

The nextflow workflow script can be ran by using the following command (activate pbs if running on the UTS HPC)
```
assembly_map_vcall.nf -c assembly_map_vcall.config \ 
    --shortreads path/to/shortreadstable.csv \
    --longreads path/to/longreadstable.csv \
    --contig path/to/assembly_edited.fa --contigname locA \
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
    * flye=2.7
    * ntedit=1.3.2