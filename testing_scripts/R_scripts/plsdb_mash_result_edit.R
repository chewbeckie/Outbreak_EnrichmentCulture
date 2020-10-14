#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)

if (length(args)==0) {
  stop("Three argument must be supplied: mash dist result, reference info file, and output file name", call.=FALSE)
}
cat("mash dist result file=",args[1], "\n",
    "reference PlsDB info file=",args[2], "\n",
    "output name=",args[3], "\n")

if(!require('tidyverse')){
  install.packages('tidyverse')
  library(tidyverse)
}

mash_output <- args[1] #mash results using plsdb.msh from PLSDB as reference
pls_info <- args[2] #plsdb.tsv downloaded from PLSDB
output_name <- args[3] 

mash <- read.table(mash_output)
colnames(mash) <- c("ref", "query", "dist", "pvalue", "sharedhashes")

pls <- read_tsv(pls_info) %>% 
  rename("ref" = ACC_NUCCORE)
              
mash.filter <- mash %>% 
  filter(pvalue < 0.1) %>%
  filter(dist < 0.1) %>% 
  right_join(., pls)

write_tsv(mash.filter, output_name, col_names = T)

cat("writing results to", output_name, "\n")