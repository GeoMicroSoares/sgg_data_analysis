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
#THIS ATE UP ~99GB RAM, NOT ADDING SPECIES
# SILVA_132_taxa <- addSpecies(SILVA_132_tax, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SILVA_132/silva_species_assignment_v132.fa.gz")
# # Write to disk
saveRDS(SILVA_132_tax, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/SILVA_132_tax_final.rds")

#Check SV table size with phyloseq
#taxa_are_rows bc of makeSequenceTable()
sgg_metadata<-read.csv("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_metadata/SGG_metadata.csv", header=TRUE)
#reorder sgg_metadata based on SGG_SVt
order<-row.names(SGG_SVt)
sgg_metadata_ord<-sgg_metadata[match(order, sgg_metadata$Sample_Code),]
rownames(sgg_metadata_ord) <- sgg_metadata_ord$Sample_Code
ps_132 <- phyloseq(otu_table(SGG_SVt, taxa_are_rows=FALSE),
               sample_data(sgg_metadata_ord),
               tax_table(SILVA_132_tax))
ps_132

ps_132.a = subset_samples(ps_132, Site.name != "Taffs Well PUMPED")
ps_132.b = subset_samples(ps_132.a, Sample_type != "Control")
ps_132.b <- prune_taxa(taxa_sums(ps_132.b) > 0, ps_132.b)
#order months
sample_data(ps_132.b)$Month = factor(sample_data(ps_132.b)$Month, 
                                 levels = c("April","August","December"))
#to relative abundances
ps_132.b = subset_taxa(ps_132.b, Kingdom %in% c("Archaea", "Bacteria"))
ps_132.b.r <-  transform_sample_counts(ps_132.b, function(x) {x/sum(x)} ) 

#Number of SVs per phyl.132um
sv.phyl.132<-as.data.frame(table(tax_table(ps_132.b.r)[, "phyl.132um"], exclude = NULL))
sv.phyl.132.ord <- sv.phyl.132[order(-sv.phyl.132$Freq),] 
sv.phyl.132.ord