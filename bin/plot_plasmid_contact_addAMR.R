#!/usr/bin/env Rscript
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

#input result from hic plasmid_map.py
asdf<-read.table("tax_links.tsv",head=F,sep='\t')
asdf$V3[is.infinite(asdf$V3)]<-max(asdf$V3[!is.infinite(asdf$V3)])

#process hic plasmid map result into matrix for making heatmap
asdf_edit <- asdf %>% 
  mutate(neg_log10_pvalue = log10(exp(V3))) %>% 
  select(2,6,neg_log10_pvalue) %>% 
  spread(key = V2, value = neg_log10_pvalue, fill = 0.0)

asdf_mt <- as.matrix(asdf_edit[,-1])
rownames(asdf_mt) <- asdf_edit$V6

#------------------------------------------------------------------------------
#input amr result from abricate
amr<-read.table("abricate_resfinder_SH.tsv",sep='\t', fill = T) %>% 
  select(V2, V14, V15) %>% 
  rename("V6" = V2)

#input plsdb mash result
pls<-read_tsv("Signf_PlsDBMashDist_SH_Oct2020.tsv", col_names = T) %>% 
  select(ref, pmlst, query)

#link amr and plasmid results together
key<-select(asdf, V1, V6) %>% 
  unique() %>% 
  left_join(.,amr, by = "V6") %>%
  set_names(c("plasmid", "contig", "AMR", "resistance")) %>% 
  group_by(contig) %>% 
  summarise(plasmid = plasmid,
            AMR_genes = paste(AMR, collapse = ", "),
            resistances = paste(resistance, collapse = ", ")) %>%
  na_if("NA") %>%  na_if("") %>% unique()

#rename some of the uncomprehensible columns
key$resistances[key$contig == "edge_555"] <- "Multi"

write.csv(key, "hicheatmap_anno_key.csv", row.names = F)

#***added types of plasmids (by web search against Pubmlst pMLST)***
#***will need to modify this part to make typing automatic instead***
key <- read.csv("hicheatmap_anno_key+pmlst.csv")

#----------------------------------------------------------------------------
#define heatmap color
col_fun <- colorRamp2(c(1:9), brewer.pal(9, "PuBu"))

#define amr annotation color
amr_col=c(brewer.pal(4,"Set1"))
names(amr_col) <- unique(key$AMR_genes)[-1]

#define resistance annotation color
res_col=c(brewer.pal(3, "Set1"))
names(res_col) <- unique(key$resistances)[-1]

#define pls type annotation color
pls_col=c(brewer.pal(8, "Pastel1"))
names(pls_col) <-sort(unique(key$type))

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
hp <- Heatmap(t(asdf_mt), name = "-log10(pvalue)", 
        show_row_dend = F, show_column_dend = F, 
        col = col_fun,
        row_names_side = "left",
        bottom_annotation = anno,
        column_labels = column_labels, column_names_rot = -45,
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(direction = "horizontal"))

#draw heatmap with legends
ComplexHeatmap::draw(hp, merge_legend = TRUE,
                     heatmap_legend_side = "right",
                     annotation_legend_side = "right")
