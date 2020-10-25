#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

if (length(args)==0) {
  stop("Three argument must be supplied: mash dist result, reference info file, and output file name", call.=FALSE)
}
cat("mash dist result file=",args[1], "\n",
    "reference chromosome genome info file=",args[2], "\n",
    "output name=",args[3], "\n")

if(!require('tidyverse')){
  install.packages('tidyverse')
  library(tidyverse)
}

mash_output <- args[1] #mash dist result "tmp_results/BacChrMashDist_SH_Oct2020.tsv"
chr_info <- args[2] #complete genome info corresponding to msh "complete_prokaryotes_enterobacteriaceae.csv"
output_name <- args[3] #output name "tmp_results/signf_BacChrMashDist_SH_Oct2020.tsv"

mash <- read.table(mash_output)
colnames(mash) <- c("ref", "query", "dist", "pvalue", "sharedhashes")

chr <- read.csv(chr_info)

chr.filter <- chr %>% 
  filter(str_detect(Replicons, 'chromosome')) %>% 
  separate(Replicons, c("Replicons_name", "genbank_refseq_accession"), sep = ":") %>% 
  mutate_at(vars(genbank_refseq_accession), ~sub('/.*', '', .)) %>% 
  mutate_at(vars(genbank_refseq_accession), ~sub(';.*', '', .)) %>% 
  mutate(query = paste(genbank_refseq_accession,".fa", sep = "")) %>% 
  select(Organism.Name, Strain, Assembly, Level, query) %>% 
  glimpse()
              
mash.filter <- mash %>% 
  separate(ref, c("remove", "ref"), sep="/") %>% 
  select(-remove) %>% 
  filter(pvalue < 0.1) %>%
  filter(dist < 0.1) %>% 
  right_join(., chr.filter)

write_tsv(mash.filter, output_name, col_names = T)

cat("writing results to", output_name, "\n")