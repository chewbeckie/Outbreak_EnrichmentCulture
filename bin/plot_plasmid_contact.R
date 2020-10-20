#!/usr/bin/env Rscript

library(ggplot2)
asdf<-read.table("tax_links.tsv",head=F,sep='\t')
asdf$V3[is.infinite(asdf$V3)]<-max(asdf$V3[!is.infinite(asdf$V3)])
neg_log10_pvalue<-log10(exp(asdf$V3))
neg_log10_pvalue[neg_log10_pvalue>7.5]<-7.5
neg_log10_pvalue[neg_log10_pvalue<=1]<-0 # these are so insignificant that it may be better to avoid visual clutter / confusion?
#asdf$V3[asdf$V3==0]=1
#contacts<-asdf$V3
(p <- ggplot(asdf, aes(V1, V2)) + geom_tile(aes(fill = neg_log10_pvalue),
     colour = "white") + scale_fill_gradient(low = "white",
     high = "steelblue"))
# "#e0e0ff",
base_size <- 14
p + theme_grey(base_size = base_size) + labs(x = "",
     y = "") + scale_x_discrete(expand = c(0, 0)) +
     scale_y_discrete(expand = c(0, 0)) + theme(
     axis.text.x = element_text(size = base_size *
         0.8, angle = 330, hjust = 0, colour = "grey50"),
     axis.text.y = element_text(size = base_size *
         1.4))


