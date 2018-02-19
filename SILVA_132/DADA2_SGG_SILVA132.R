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
#### quick transformations ####
ps_132.a = subset_samples(ps_132, Site.name != "Taffs Well PUMPED")
ps_132.b = subset_samples(ps_132.a, Sample_type != "Control")
ps_132.b <- prune_taxa(taxa_sums(ps_132.b) > 0, ps_132.b)
#order months
sample_data(ps_132.b)$Month = factor(sample_data(ps_132.b)$Month, 
                                 levels = c("April","August","December"))
#to relative abundances
ps_132.b = subset_taxa(ps_132.b, Kingdom %in% c("Archaea", "Bacteria"))
ps_132.b.r <-  transform_sample_counts(ps_132.b, function(x) {x/sum(x)} ) 

#Number of SVs per phylum ####
sv.phyl.132<-as.data.frame(table(tax_table(ps_132.b.r)[, "Phylum"], exclude = NULL))
sv.phyl.132.ord <- sv.phyl.132[order(-sv.phyl.132$Freq),] 
sv.phyl.132.ord
head(sv.phyl.132.ord)
#unclassified SV only go down by 0.4% with SILVA_132
#first 4 classified phyla's SVs make up for 59.06% of the dataset.

#### taxonomy patterns ####
#Check how tax is doing
ps_132.b.r.glom <- tax_glom(ps_132.b.r, taxrank = 'Phylum', NArm = FALSE)
# create dataframe from phyloseq object
ps_132.b.r.glom.ps_132df <- data.table(psmelt(ps_132.b.r.glom))
# convert Phylum to a character vector from a factor because R
ps_132.b.r.glom.ps_132df$Phylum <- as.character(ps_132.b.r.glom.ps_132df$Phylum)
# group dataframe by Phylum, calculate median rel. abundance
ps_132.b.r.glom.ps_132df[, median := median(Abundance, na.rm = TRUE), 
                 by = "Phylum"]
ps_132.b.r.glom.ps_132df.noProteo<-ps_132.b.r.glom.ps_132df[!grepl("Proteobacteria", ps_132.b.r.glom.ps_132df$Phylum),]
#New table, with Proteobacterial classes only
ps_132.b.r.ProtCl = subset_taxa(ps_132.b.r, Phylum == "Proteobacteria")
ps_132.b.r.ProtCl.glom <- tax_glom(ps_132.b.r.ProtCl, taxrank = 'Class', NArm = FALSE)
# create dataframe from phyloseq object
ps_132.b.r.ProtCl.glom.ps_132df <- data.table(psmelt(ps_132.b.r.ProtCl.glom))
# convert Class to a character vector from a factor because R
ps_132.b.r.ProtCl.glom.ps_132df$Class <- as.character(ps_132.b.r.ProtCl.glom.ps_132df$Class)
# group dataframe by Phylum, calculate median rel. abundance
ps_132.b.r.ProtCl.glom.ps_132df[, median := median(Abundance, na.rm = TRUE), 
                        by = "Class"]
# Change name to remainder of Phylum less than 1%
ps_132.b.r.ProtCl.glom.ps_132df[(median <= 0.01), Class := "Other Proteobacteria"]
#join the two data frames
ps_132.b.r.ProtCl.glom.ps_132df.noP <- subset(ps_132.b.r.ProtCl.glom.ps_132df, select = -Phylum)
#same name= allows merging
names(ps_132.b.r.glom.ps_132df.noProteo)[names(ps_132.b.r.glom.ps_132df.noProteo) == 'Phylum'] <- 'Taxa'
names(ps_132.b.r.ProtCl.glom.ps_132df.noP)[names(ps_132.b.r.ProtCl.glom.ps_132df.noP) == 'Class'] <- 'Taxa'
#join the two
ps_132.b.r.ProteoCl.AllPhy <- rbind(ps_132.b.r.glom.ps_132df.noProteo, ps_132.b.r.ProtCl.glom.ps_132df.noP)
taxa_number<-length(unique(ps_132.b.r.ProteoCl.AllPhy$Taxa))
ps_132.b.r.ProteoCl.AllPhy.plot<-ggplot(ps_132.b.r.ProteoCl.AllPhy[Abundance > 0], 
                                    aes(x = Site.name, y = Abundance/3, fill = Taxa)) + 
  facet_grid(Month~.) +
  geom_bar(stat = "identity") +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(nrow = taxa_number/2)) +
  ylab("Relative Abundance\n") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(taxa_number), 
                    na.value = "gray48")
#totals vary because the bars = medians
ps_132.b.r.ProteoCl.AllPhy.plot

#### going deep ####
ps.b.deep_132 = subset_taxa(ps_132.b.r, 
                            Phylum=="Epsilonbacteraeota" | Phylum=="Proteobacteria" | Phylum=="Patescibacteria" | Phylum=="Omnitrophicaeota")
ps.b.deep_132.glom <- tax_glom(ps.b.deep_132, taxrank = 'Genus', NArm = FALSE)
ps.b.deep_132.glom.df <- data.table(psmelt(ps.b.deep_132.glom))
colourCount = length(unique(ps.b.deep_132.glom.df[Abundance > 0.03]$Genus))+1
ps.b.deep_132.glom.df.plot<-ggplot(ps.b.deep_132.glom.df[Abundance > 0.03], 
                                        aes(x = Site.name, y = Abundance/3, fill = Genus)) + 
  facet_grid(Month~.) +
  geom_bar(stat = "identity") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 1)) + 
  guides(fill = guide_legend(nrow = colourCount)) +
  ylab("Relative Abundance\n") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(colourCount), 
                    na.value = "gray48")
ps.b.deep_132.glom.df.plot
