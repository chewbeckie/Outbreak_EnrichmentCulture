#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

#------------------------------------------------------------------------------
#define input and output here
setwd("contigs_annotate_results_20201025/locB_results/")

#hic plasmid map.py result
asdf_path <- "~/github/Outbreak_EnrichmentCulture/data/results/tax_links_locB.tsv"
#ABRicate AMR gene screening result
amr_path <- "locB_resfinder_amr.tsv"
#ABRicate pMLST gene screening result
pmlst_path <- "locB_plsfinder_pmlst.tsv"
#output summary path (as tsv)
output_path <- "locB_amr_pmlst_plasmid_summary.tsv"

#------------------------------------------------------------------------------
#input result from hic plasmid_map.py
asdf<-read.table(asdf_path,head=F,sep='\t')
asdf$V3[is.infinite(asdf$V3)]<-max(asdf$V3[!is.infinite(asdf$V3)])

#process hic plasmid map result into matrix for making heatmap
asdf <- asdf %>% 
  mutate(neg_log10_pvalue = log10(exp(V3))) %>% 
  #filter pvalue for visualization
  mutate(neg_log10_pvalue = case_when(neg_log10_pvalue > 7.5 ~ 7.5,
                                      neg_log10_pvalue <= 1  ~ 0,
                                      TRUE ~ neg_log10_pvalue)) %>% 
  filter(neg_log10_pvalue > 0) 
  
asdf_edit <- asdf %>% 
  select(2,6,neg_log10_pvalue) %>%
  set_names(c("genus", "contig", "neg_log10_pvalue")) %>% 
  spread(key = genus, value = neg_log10_pvalue, fill = 0.0)

asdf_mt <- as.matrix(asdf_edit[,-1])
rownames(asdf_mt) <- asdf_edit$contig

#------------------------------------------------------------------------------
#input amr result from abricate
amr<-read.table(amr_path,sep='\t', fill = T) %>% 
  set_names(c("FILE","SEQUENCE","START","END","STRAND","GENE","COVERAGE",
              "COVERAGE_MAP","GAPS","percentCOVERAGE","percentIDENTITY",
              "DATABASE","ACCESSION","AMR","RESISTANCE")) %>% 
  filter(percentIDENTITY >= 90 & percentCOVERAGE >=90) %>% 
  select(SEQUENCE, AMR, RESISTANCE) %>% 
  rename("contig" = SEQUENCE)

#input plsfinder result from abricate
plsfinder<- read.table(pmlst_path,sep='\t', fill = T) %>% 
  set_names(c("FILE","SEQUENCE","START","END","STRAND","pMLST","COVERAGE",
              "COVERAGE_MAP","GAPS","percentCOVERAGE","percentIDENTITY",
              "DATABASE","ACCESSION","PRODUCT","RESISTANCE")) %>% 
  filter(percentIDENTITY >= 90 & percentCOVERAGE >=90) %>% 
  select(SEQUENCE, pMLST) %>% 
  rename("contig" = SEQUENCE)

#link amr, plsfinder and plasmid results together
key<-select(asdf, V1, V6) %>% 
  set_names(c("plasmid", "contig")) %>% 
  unique() %>% 
  left_join(.,amr, by = "contig") %>%
  left_join(.,plsfinder, by = "contig") %>% 
  group_by(contig) %>% 
  summarise(plasmid = plasmid,
            AMR_genes = paste(AMR, collapse = ";"),
            resistances = paste(RESISTANCE, collapse = ";"),
            type = paste(pMLST)) %>%
  mutate_at(vars(AMR_genes, resistances), ~sub('.*NA','NA',.)) %>% 
  na_if(";") %>% na_if("NA") %>% 
  unique()

#write annotation to output
write.table(key, output_path, row.names = F, sep = "\t")

#additional------------------------------------------------------------------
#rename some of the uncomprehensible columns
#key$resistances[key$contig == "edge_555"] <- "Multi"

#----------------------------------------------------------------------------
#define heatmap color
col_fun <- colorRamp2(c(1:7), brewer.pal(7, "PuBu"))

#define amr annotation color
ncol<-length(unique(key$AMR_genes)[-1])
amr_col<-c(brewer.pal(ncol,"Set1"))
names(amr_col) <- unique(key$AMR_genes)[-1]

#define resistance annotation color
ncol<-length(unique(key$resistances)[-1])
res_col<-c(brewer.pal(ncol, "Set1"))
names(res_col) <- unique(key$resistances)[-1]

#define pls type annotation color
ncol<-length(unique(key$type)[-1])
pls_col=c(brewer.pal(ncol, "Set3"))
names(pls_col) <-sort(unique(key$type)[-1])

#define annotation
anno <- HeatmapAnnotation(AMRgene = key$AMR_genes,
                          Antibiotics = key$resistances,
                          Type = key$type,
                          col = list(AMRgene = amr_col,
                                     Antibiotics = res_col,
                                     Type = pls_col),
                          na_col = "lightgrey")

#use plasmid name as column names
column_labels <- structure(key$plasmid, names = key$contig)

#draw heatmap
Heatmap(t(asdf_mt), name = "-log10(pvalue)", 
        show_row_dend = F, show_column_dend = F, 
        col = col_fun,
        row_names_side = "left",
        bottom_annotation = anno,
        column_labels = column_labels, column_names_rot = -45,
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(direction = "vertical"))