# assembly-metaflye_and_polish.nf usage

De novo metagenome assembly using ONT long-read data and correction with illumina short-read data.

1. Prepare a `.csv` runtable with information about the sample name, path to read1 and path to read 2 of each individual Illumina shortread dataset

```
sampleId,read1,read2
sample_1,path/to/r1,path/to/r2
sample_2,anotherpath/to/r1,anotherpath/to/r2
```

2. Similarly, prepare a `.csv` runtable with information about the sample name and path to the ONT longread dataset

```
sampleId,longread
sample_a,path/to/longread
```

3. specify the genome size (`params.gs`) in the `assembly-metaflye_and_polish.nf`, the genome size is set at 100m for the generation of WG/SH assemblies.

4. *(optional)* specify directory locations for input folder and output folder in `assembly-metaflye_and_polish.config`.

5. run the nextflow workflow script by using the following command (activate pbs if running on the UTS HPC).
```
assembly-metaflye_and_polish.nf -c assembly-metaflye_and_polish.config \
 --shortreads shortreads.csv --longreads longreads.csv \
 --contigname SH \
 -profile conda,pbs
```

## Dependence
* nextflow
* conda (dependencies are listed in `assembly-metaflye_and_polish.config`)
    * flye=2.7
    * fastp
    * ntedit=1.3.2
* python