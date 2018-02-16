#completely adapted from https://benjjneb.github.io/dada2/tutorial.html
library(dada2); packageVersion("dada2")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")
library(data.table); packageVersion("data.table")
library(DESeq2); packageVersion("DESeq2")
library(dplyr)
library(plyr)
library(magrittr)
library(RColorBrewer)
library(phangorn)
library(DECIPHER)
library("structSSI")
library(ips)
library(msa)
library(picante)
library(igraph)
library(cowplot)

#first check that seqtab_final corresponds to version corresponding to original QC parameters
SGG_SVt = readRDS("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/seqtab_final.rds")
SILVA_132_tax <- assignTaxonomy(SGG_SVt, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SILVA_132/silva_nr_v132_train_set.fa.gz", 
                      multithread=TRUE)
SILVA_132_taxa <- addSpecies(SILVA_132_tax, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SILVA_132/silva_species_assignment_v132.fa.gz")
# Write to disk
saveRDS(SILVA_132_taxa, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/SILVA_132_taxa_final.rds")

#Check SV table size with phyloseq
#taxa_are_rows bc of makeSequenceTable()
sgg_metadata<-read.csv("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_metadata/SGG_metadata.csv", header=TRUE)
#reorder sgg_metadata based on seqtab_bimR
order<-row.names(seqtab_bimR)
sgg_metadata_ord<-sgg_metadata[match(order, sgg_metadata$Sample_Code),]
rownames(sgg_metadata_ord) <- sgg_metadata_ord$Sample_Code
ps_132 <- phyloseq(otu_table(seqtab_bimR, taxa_are_rows=FALSE),
               sample_data(sgg_metadata_ord),
               tax_table(SILVA_132_taxa))
ps_132