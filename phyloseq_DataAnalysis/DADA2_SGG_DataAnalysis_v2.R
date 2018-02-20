# Load libraries ----
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
library(plotly)
library(ggnet)
library(network)
library(sna)
library(intergraph)

# Import SV table, taxonomy ----
SGG_SVt = readRDS("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/seqtab_final.rds")
SGG_Tax = readRDS("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/tax_final.rds")
#Check SV table size with phyloseq
#taxa_are_rows bc of makeSequenceTable()
sgg_metadata<-read.csv("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_metadata/SGG_metadata.csv", header=TRUE)
#reorder sgg_metadata based on seqtab_bimR
order<-row.names(SGG_SVt)
sgg_metadata_ord<-sgg_metadata[match(order, sgg_metadata$Sample_Code),]
rownames(sgg_metadata_ord) <- sgg_metadata_ord$Sample_Code

ps <- phyloseq(otu_table(SGG_SVt, taxa_are_rows=FALSE),
               sample_data(sgg_metadata_ord), 
               tax_table(SGG_Tax))
#this is the full dataset
ps

# 0 - Basic stats, QC ----
#how many reads left in seqtab_bimR?
sum(rowSums(otu_table(ps)))
#How many OTUs in ps?
ntaxa(otu_table(ps))
#70,937 SVs to 9,656,492 reads

#Seq Depth
psdt = data.table(as(sample_data(ps), "data.frame"),
                  TotalReads = sample_sums(ps), keep.rownames = TRUE,
                  Batch = sgg_metadata_ord$Batch)
setnames(psdt, "rn", "SampleID")
pSeqDepth = ggplot(psdt, aes(TotalReads)) + 
  geom_histogram(binwidth = 5000) + 
  ggtitle("Sequencing Depth") +
  facet_grid(.~Batch)
pSeqDepth
ggsave("pSeqDepth.png", width=14, height=6)

##0.5 - QC, contaminant control ====
readsumsdf_samples<-data.frame(nreads = sample_sums(ps),
                               samples = rownames(as.data.frame(sample_sums(ps))),
                               batch = sample_data(ps)$Batch)
readsumsdf_samples_ord<-readsumsdf_samples[order(readsumsdf_samples$samples),]

samples_table <- table(as.numeric(readsumsdf_samples$samples))
samples_levels <- names(samples_table)[order(samples_table)]
samples_levels[c(121,122,123,124,125,126,127,128)] <- c("P1","P2","P3","P4","S1","S2","S3","S4")
readsumsdf_samples$samples2 <- factor(readsumsdf_samples$samples, levels = samples_levels)

op=ggplot(readsumsdf_samples, aes(x=samples2,y=nreads, 
                                  color=batch, fill=batch)) +
  geom_bar(stat = "identity") +
  scale_y_log10(breaks=c(0,100,500,1000,5000,10000,25000,50000,100000,200000)) +
  geom_hline(yintercept =1416,
             color="dark red",
             linetype="dashed") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
op
ggsave("SGG_reads_batches.png", width=20, height=10)

##Analyse SVs in controls
ps.control = subset_samples(ps, Sample_type == "Control")
ps.control.f = prune_taxa(taxa_sums(ps.control) > 0, ps.control)
plot_bar(ps.control.f, fill="Genus") +
  annotate("text", y= 150, x ="P4",label=paste("Total SVs: 12\nTotal Reads:742")
           ,hjust=0.5, size=6)
ggsave("SGG_controlSVs.png", width=12, height=6)

control.taxa = names(taxa_sums(ps.control.f))
ps.b.r.control.taxa = prune_taxa(control.taxa, ps.r)
ps.b.r.control.glom <- tax_glom(ps.b.r.control.taxa, taxrank = 'Genus', NArm = FALSE)
ps.b.r.control.glom.psdf <- data.table(psmelt(ps.b.r.control.glom))
ps.b.r.control.glom.psdf$Month = factor(ps.b.r.control.glom.psdf$Month, 
                                        levels = c("April","August","December","Control"))
ps.b.r.control.glom.psdf$Genus <- as.character(ps.b.r.control.glom.psdf$Genus)
ps.b.r.control.glom.psdf[, median := median(Abundance, na.rm = TRUE), 
                         by = "Genus"]
ps.b.r.control.glom.psdf.plot<-ggplot(ps.b.r.control.glom.psdf[Abundance > 0], 
                                      aes(x = Site.name, y = Abundance/3, fill = Genus)) + 
  facet_grid(Month~.) +
  # scale_y_log10(breaks=c(10,100,1000,30000)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance\n")
ps.b.r.control.glom.psdf.plot
ggsave("SGG_controlSVs_overall_relabund.png", width=20, height=10)
#absolute abundances
ps.control.taxa = prune_taxa(control.taxa, ps)
ps.control.glom <- tax_glom(ps.control.taxa, taxrank = 'Genus', NArm = FALSE)
ps.control.glom.psdf <- data.table(psmelt(ps.control.glom))
ps.control.glom.psdf$Month = factor(ps.control.glom.psdf$Month, 
                                    levels = c("April","August","December","Control"))
ps.control.glom.psdf$Genus <- as.character(ps.control.glom.psdf$Genus)
ps.control.glom.psdf[, median := median(Abundance, nam = TRUE), 
                     by = "Genus"]
ps.control.glom.psdf.plot<-ggplot(ps.control.glom.psdf[Abundance > 0], 
                                  aes(x = Site.name, y = Abundance/3, fill = Genus)) + 
  facet_grid(Month~.) +
  # scale_y_log10(breaks=c(10,100,1000,30000)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance\n")
ps.control.glom.psdf.plot
ggsave("SGG_controlSVs_overall_absabund.png", width=20, height=10)

#0.5 - Custom datasets  ----
#now without Controls and Taffs Well PUMPED
ps.a = subset_samples(ps, Site.name != "Taffs Well PUMPED")
ps.b = subset_samples(ps.a, Sample_type != "Control")
ps.b <- prune_taxa(taxa_sums(ps.b) > 0, ps.b)
#order months
sample_data(ps.b)$Month = factor(sample_data(ps.b)$Month, 
                                 levels = c("April","August","December"))
#to relative abundances
ps.b = subset_taxa(ps.b, Kingdom %in% c("Archaea", "Bacteria"))
ps.b.r <-  transform_sample_counts(ps.b, function(x) {x/sum(x)} ) 

#Number of SVs per phylum
sv.phyl<-as.data.frame(table(tax_table(ps.b.r)[, "Phylum"], exclude = NULL))
sv.phyl.ord <- sv.phyl[order(-sv.phyl$Freq),] 
sv.phyl.ord
#Number of sequences per SV, Sequences per sample
readsumsdf = data.frame(nreads = sort(taxa_sums(ps.b), TRUE), sorted = 1:ntaxa(ps.b), type = "SVs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps.b), 
                                                        TRUE), sorted = 1:nsamples(ps.b), type = "Samples"))
p_read_sums = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p_read_sums + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

# 1 - Taxonomic analysis ----
# Phylum-level tax, transformed ====
# agglomerate taxa
ps.b.r.glom <- tax_glom(ps.b.r, taxrank = 'Phylum', NArm = FALSE)
# create dataframe from phyloseq object
ps.b.r.glom.psdf <- data.table(psmelt(ps.b.r.glom))
# convert Phylum to a character vector from a factor because R
ps.b.r.glom.psdf$Phylum <- as.character(ps.b.r.glom.psdf$Phylum)
# group dataframe by Phylum, calculate median rel. abundance
ps.b.r.glom.psdf[, median := median(Abundance, na.rm = TRUE), 
                 by = "Phylum"]
# Change name to remainder of Phylum less than 1%
ps.b.r.glom.psdf[(median <= 0.01), Phylum := "Other"]

ps.b.r.glom.psdf.plot<-ggplot(ps.b.r.glom.psdf[Abundance > 0], 
                              aes(x = Site.name, y = Abundance/3, fill = Phylum)) + 
  facet_grid(Month~.) +
  geom_bar(stat = "identity") +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  # coord_cartesian(ylim = c(0, 1)) +
  ylab("Relative Abundance\n")
#totals vary because the bars = medians
ps.b.r.glom.psdf.plot
ggsave("SGG-0.01pc-Phylum-R-Month+Site.plot.png", height=8, width=16)

# Phylum-, Proteobacterial class-level tax, not transformed ====
#Removing Proteobacteria from previous table
ps.b.r.glom.psdf.noProteo<-ps.b.r.glom.psdf[!grepl("Proteobacteria", ps.b.r.glom.psdf$Phylum),]
#New table, with Proteobacterial classes only
ps.b.r.ProtCl = subset_taxa(ps.b.r, Phylum == "Proteobacteria")
ps.b.r.ProtCl.glom <- tax_glom(ps.b.r.ProtCl, taxrank = 'Class', NArm = FALSE)
# create dataframe from phyloseq object
ps.b.r.ProtCl.glom.psdf <- data.table(psmelt(ps.b.r.ProtCl.glom))
# convert Class to a character vector from a factor because R
ps.b.r.ProtCl.glom.psdf$Class <- as.character(ps.b.r.ProtCl.glom.psdf$Class)
# group dataframe by Phylum, calculate median rel. abundance
ps.b.r.ProtCl.glom.psdf[, median := median(Abundance, na.rm = TRUE), 
                        by = "Class"]
# Change name to remainder of Phylum less than 1%
ps.b.r.ProtCl.glom.psdf[(median <= 0.01), Class := "Other Proteobacteria"]
#Check if they're all there
# ggplot(ps.b.r.ProtCl.glom.psdf[Abundance > 0], 
#        aes(x = Site.name, y = Abundance/3, fill = Class)) +
#   facet_grid(Month~.) +
#   geom_bar(stat = "identity")
#join the two data frames
ps.b.r.ProtCl.glom.psdf.noP <- subset(ps.b.r.ProtCl.glom.psdf, select = -Phylum)
#same name= allows merging
names(ps.b.r.glom.psdf.noProteo)[names(ps.b.r.glom.psdf.noProteo) == 'Phylum'] <- 'Taxa'
names(ps.b.r.ProtCl.glom.psdf.noP)[names(ps.b.r.ProtCl.glom.psdf.noP) == 'Class'] <- 'Taxa'
#join the two
ps.b.r.ProteoCl.AllPhy <- rbind(ps.b.r.glom.psdf.noProteo, ps.b.r.ProtCl.glom.psdf.noP)

ps.b.r.ProteoCl.AllPhy.plot<-ggplot(ps.b.r.ProteoCl.AllPhy[Abundance > 0], 
                                    aes(x = Site.name, y = Abundance/3, fill = Taxa)) + 
  facet_grid(Month~.) +
  geom_bar(stat = "identity") +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  # coord_cartesian(ylim = c(0, 1)) +
  ylab("Relative Abundance\n")
#totals vary because the bars = medians
ps.b.r.ProteoCl.AllPhy.plot
ggsave("SGG-0.01pc-Phylum+ProtoCl-R-Month+Site.plot.png", height=8, width=16)

# Genus-level tax, transformed ====
#Tracking taxa responsible for Beta, Epsilonproteobacteria dominances
# Gallionellaceae, not transformed ####
ps.b.gallionelac = subset_taxa(ps.b, Family=="Gallionellaceae")
ps.b.gallionelac.plot<-plot_bar(ps.b.gallionelac, 
                                fill="Genus", facet_grid=~Month, x="Site.name") + 
  geom_bar(aes(color=Genus, fill=Genus), 
           stat="identity", position="stack")
ps.b.gallionelac.plot

ps.b.r.galliGen = subset_taxa(ps.b.r, Genus=="Gallionella" | Genus=="Sideroxydans" | Genus=="Ferriphaselus")
ps.b.r.galliGen.glom <- tax_glom(ps.b.r.galliGen, taxrank = 'Genus', NArm = FALSE)
ps.b.r.galliGen.glom.psdf <- data.table(psmelt(ps.b.r.galliGen.glom))
ps.b.r.galliGen.glom.psdf$Genus <- as.character(ps.b.r.galliGen.glom.psdf$Genus)

ps.b.r.galliGen.glom.psdf.plot<-ggplot(ps.b.r.galliGen.glom.psdf,
                                       aes(x=Site.name,y=Abundance,fill=Genus)) + 
  geom_boxplot() + 
  facet_grid(Month~.)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))
ps.b.r.galliGen.glom.psdf.plot
ggsave("SGG-GallionellaceaeGenera-R-Month-Site.plot.png",
       width = 14,
       height = 6)
# Campylobacterales, not transformed ####
ps.b.helicob = subset_taxa(ps.b, Family=="Helicobacteraceae")
ps.b.helicob.plot<-plot_bar(ps.b.helicob, 
                            fill="Genus", facet_grid=~Month, x="Site.name") + 
  geom_bar(aes(color=Genus, fill=Genus), 
           stat="identity", position="stack")
ps.b.helicob.plot

ps.b.r.HelicobGen = subset_taxa(ps.b.r, Genus=="Sulfuricurvum" | Genus=="Sulfurimonas" | Genus=="Sulfurovum")
ps.b.r.HelicobGen.glom <- tax_glom(ps.b.r.HelicobGen, taxrank = 'Genus', NArm = FALSE)
ps.b.r.HelicobGen.glom.psdf <- data.table(psmelt(ps.b.r.HelicobGen.glom))
ps.b.r.HelicobGen.glom.psdf$Genus <- as.character(ps.b.r.HelicobGen.glom.psdf$Genus)

ps.b.r.HelicobGen.glom.psdf.plot<-ggplot(ps.b.r.HelicobGen.glom.psdf,
                                         aes(x=Site.name,y=Abundance,fill=Genus)) + 
  geom_boxplot() + 
  facet_grid(Month~.)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))
ps.b.r.HelicobGen.glom.psdf.plot
ggsave("SGG-HelicobacteraceaeGenera-R-Month-Site.plot.png",
       width = 14,
       height = 6)

#How about all these together?
ps.b.r.HelicobGall = subset_taxa(ps.b.r, Genus=="Sulfuricurvum" | Genus=="Sulfurimonas" | Genus=="Sulfurovum" | Genus=="Gallionella" | Genus=="Sideroxydans" | Genus=="Ferriphaselus")
ps.b.r.HelicobGall.glom <- tax_glom(ps.b.r.HelicobGall, taxrank = 'Genus', NArm = FALSE)
ps.b.r.HelicobGall.glom.psdf <- data.table(psmelt(ps.b.r.HelicobGall.glom))
ps.b.r.HelicobGall.glom.psdf <- data.frame(lapply(ps.b.r.HelicobGall.glom.psdf, function(x) {gsub("Crumlin Navigation", "Crum. Nav.", x)}))
ps.b.r.HelicobGall.glom.psdf$Genus <- as.character(ps.b.r.HelicobGall.glom.psdf$Genus)
ps.b.r.HelicobGall.glom.psdf$Abundance <- as.numeric(as.character(ps.b.r.HelicobGall.glom.psdf$Abundance))

ps.b.r.HelicobGall.glom.psdf.plot<-ggplot(ps.b.r.HelicobGall.glom.psdf,
                                          aes(x=Genus,y=Abundance,fill=Genus)) + 
  geom_point(aes(color=Genus, fill=Genus,
                 size=Abundance)) + 
  facet_grid(Month~Site.name)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ps.b.r.HelicobGall.glom.psdf.plot
ggsave("SGG-Helicobacteraceae+Gallionellaceae-Genera-R-Month-Site.plot.png",
       width = 18,
       height = 8)

# Proteobacteria-less tax, transformed ====
ps.b.r.noProt<-subset_taxa(ps.b.r, Phylum!="Proteobacteria")
#check phyla
sv.phyl.noProt<-as.data.frame(table(tax_table(ps.b.r.noProt)[, "Phylum"], exclude = NULL))
sv.phyl.noProt.ord <- sv.phyl.noProt[order(-sv.phyl.noProt$Freq),] 
#let's go for Phyla w/more than 1000 SV's following Proteobacteria
ps.b.r.noProt.glom <- tax_glom(ps.b.r.noProt, taxrank = 'Phylum', NArm = FALSE)
ps.b.r.noProt.glom.psdf <- data.table(psmelt(ps.b.r.noProt.glom))
ps.b.r.noProt.glom.psdf$Phylum <- as.character(ps.b.r.noProt.glom.psdf$Phylum)
ps.b.r.noProt.glom.psdf[, median := median(Abundance, na.rm = TRUE), 
                        by = "Phylum"]
ps.b.r.noProt.glom.psdf[(median <= 0.005), Phylum := "Other"]

ps.b.r.noProt.glom.psdf.plot<-ggplot(ps.b.r.noProt.glom.psdf[Abundance > 0], 
                                     aes(x = Phylum, y = Abundance, fill = Phylum)) + 
  facet_grid(Month~Site.name) +
  geom_point(aes(color=Phylum, fill=Phylum,
                 size=Abundance)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+ 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance\n")
#totals vary because the bars = medians
ps.b.r.noProt.glom.psdf.plot
ggsave("SGG-PhylaBelowProteobacteria-R-Month+Site.plot.png", 
       height=8, 
       width=20)
#
ps.b.r.OmniParcu.class<-as.data.frame(table(tax_table(ps.b.r.OmniParcu)[, "Class"], exclude = NULL))
ps.b.r.OmniParcu.class.ord <- ps.b.r.OmniParcu.class[order(-ps.b.r.OmniParcu.class$Freq),] 
ps.b.r.OmniParcu.class.ord

ps.b.r.OmniParcu = subset_taxa(ps.b.r, 
                               Phylum=="Parcubacteria" | Phylum=="Omnitrophica" )
ps.b.r.OmniParcu.glom <- tax_glom(ps.b.r.OmniParcu, taxrank = 'Class', NArm = FALSE )
ps.b.r.OmniParcu.glom.psdf <- data.table(psmelt(ps.b.r.OmniParcu.glom))
ps.b.r.OmniParcu.glom.psdf$Class <- as.character(ps.b.r.OmniParcu.glom.psdf$Class)
# plot(sort(ps.b.r.OmniParcu.glom.psdf$Abundance, decreasing = TRUE))
#plot only classes representing >1% of the sample profile in the future
ps.b.r.OmniParcu.glom.psdf.plot<-ggplot(ps.b.r.OmniParcu.glom.psdf[Abundance > 0], 
                                        aes(x = Class, y = Abundance, fill = Class)) + 
  facet_grid(Month~Site.name, space = "free_x") +
  geom_point(aes(color=Class, fill=Class,
                 size=Abundance)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+ 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance\n")
#totals vary because the bars = medians
ps.b.r.OmniParcu.glom.psdf.plot
ggsave("SGG-PhylaBelowProteobacteria-R-Month+Site.plot.png", 
       height=8, 
       width=20)

# Lindsay, Cefn Hengoed, both have ~50% Other, transformed ====
ps.b.r.LindsayCHen = subset_samples(ps.b.r, 
                                    Site.name=="Lindsay" | Site.name=="Cefn Hengoed" )
ps.b.r.LindsayCHen.glom <- tax_glom(ps.b.r.LindsayCHen, taxrank = 'Phylum', NArm = FALSE )
ps.b.r.LindsayCHen.glom.psdf <- data.table(psmelt(ps.b.r.LindsayCHen.glom))
ps.b.r.LindsayCHen.glom.psdf$Phylum <- as.character(ps.b.r.LindsayCHen.glom.psdf$Phylum)
ps.b.r.LindsayCHen.glom.psdf[, median := median(Abundance, na.rm = TRUE), 
                             by = "Phylum"]
ps.b.r.LindsayCHen.glom.psdf[(median <= 0.0001), Phylum := "Other"]
ps.b.r.LindsayCHen.glom.psdf.plot<-ggplot(ps.b.r.LindsayCHen.glom.psdf[Abundance > 0], 
                                          aes(x = Site.name, y = Abundance/3, fill = Phylum)) + 
  facet_grid(Month~.) +
  geom_bar(stat="identity", 
           position = position_stack(),
           colour="black", size=.1) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+ 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance\n")
ps.b.r.LindsayCHen.glom.psdf.plot

ps.b.r.LindsayCHen = subset_samples(ps.b.r, Site.name=="Lindsay")
ps.b.r.glom.phy <- tax_glom(ps.b.r.LindsayCHen, taxrank = 'Phylum', NArm = FALSE )
ps.b.r.glom.phy.psdf <- data.table(psmelt(ps.b.r.glom.phy))
ps.b.r.glom.phy.psdf.f <- ps.b.r.glom.phy.psdf[ , c("Abundance","Phylum")]
ps.b.r.glom.phy.psdf.f.s<-ps.b.r.glom.phy.psdf.f[, .(count=.N, var=mean(Abundance)), by = Phylum]

ps.b.r.glom.phy.psdf.f.s <- transform(ps.b.r.glom.phy.psdf.f.s, Phylum=reorder(Phylum, -var)) 
ggplot(ps.b.r.glom.phy.psdf.f.s, aes(x=Phylum, y=var, color=Phylum)) +
  geom_point() +
  geom_hline(yintercept = 0.01,
             color="gray",
             linetype="dashed") +
  ylab="Relative abundance" +
  theme(axis.text.x = element_text(angle=40, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")

#seems about 50% of their profiles' class-level taxonomy hasn't been classified

# 3 - Ordination on ps ----
# nMDS and PCoA on dataset
# ps.b.r.dist <- phyloseq::distance(ps.b.r, method="bray")
# Calculate ordination
ps.b.r.pcoa  <- ordinate(ps.b.r, "PCoA", distance="jsd", 
                         k=2, trymax=1e3, weighted=TRUE)
#How good is this, i.e. how much of the variance is explained by the first few components?
plot_scree(ps.b.r.pcoa)
evals<-ps.b.r.pcoa$values$Eigenvalues
#Axis 1 and 2
pcoa.a12<-sum(evals[1:2])/sum(evals)
pcoa.a12.r<-round(pcoa.a12, digits=4)*100
#Axis 1 and 3
pcoa.a13<-sum(evals[c(1,3)])/sum(evals)
pcoa.a13.r<-round(pcoa.a13, digits=4)*100

ps.b.r.pcoa.plot.a1 <- plot_ordination(ps.b.r, ps.b.r.pcoa, axes = c(1,2),
                                       color="Site.name", shape="Month") +
  geom_polygon(aes(fill=Site.name), alpha=0.5) +
  geom_point(alpha=0.5, size=4)+
  scale_alpha(guide = 'none') +
  theme(legend.position = "none") +
  annotate("text", y= 0.3, x =0.25,label=paste("Axes 1 and 2:\n",pcoa.a12.r,"% variation")
           ,hjust=0.5, size=5)

ps.b.r.pcoa.plot.a2 <- plot_ordination(ps.b.r, ps.b.r.pcoa, axes = c(1,3),
                                       color="Site.name", shape="Month")+
  geom_polygon(aes(fill=Site.name), alpha=0.5) +
  geom_point(alpha=0.5, size=4) +
  scale_alpha(guide = 'none')+
  annotate("text", y= 0.2, x =0.25,label=paste("Axes 1 and 3:\n",pcoa.a13.r,"% variation")
           ,hjust=0.5, size=5)

ps.b.r.pcoa.plot<-plot_grid(ps.b.r.pcoa.plot.a1, ps.b.r.pcoa.plot.a2, 
                            labels=c("A","B"))
ps.b.r.pcoa.plot
ggsave("SGG-PCoA-Axes123-R-Month+Site.plot.png", 
       height=6, 
       width=16)

#### But how much do site and sampling month, etc affect community structure?
set.seed(1)
# Calculate bray curtis distance matrix
ps_jsd <- phyloseq::distance(ps.b.r, method = "jsd")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps.b.r))
# Adonis test: "how is variation attributed to different experimental" variables?
# null hypothesis: tested metadata categories have the same centroid.
adonis(ps_jsd ~ Site.name + Month + Batch + Temp.C + EC.uS.cm + K.µg.ml + Na.µg.ml + Cl..µg.ml + SO4..µg.ml, 
       data = sampledf, strata = sampledf$Time_point, permutations=1e3)
# Only Site name, Month, Temp.C and Cl..µg.ml explain significantly community variation
# meaning according to Site, Month Temperature and Cl individually, samples don't have a same centroid.

# Homogeneity of dispersion test for significant categories
# null hypothesis: tested metadata categories have the same dispersions (variance).
# betadisper() calculates average distances of samples to group centroids (from metadata categories)
# permutest() calculates ANOVAS on average distances from betadisper() to test if one or more groups is more variable than others
beta_s <- betadisper(ps_jsd, sampledf$Site.name)
permutest(beta_s)
beta_m <- betadisper(ps_jsd, sampledf$Month)
permutest(beta_m)
beta_t <- betadisper(ps_jsd, sampledf$Temp.C)
permutest(beta_t)
beta_c <- betadisper(ps_jsd, sampledf$Cl..µg.ml)
permutest(beta_c)
# we can reject the null hypothesis that site names have the same dispersions (variance).
# HOWEVER this might mean that adonis results can't be trusted for this category.

# Tukey HSD as an alternative to permutest()
#Create a set of confidence intervals on the differences between the means of the levels of the grouping factor
# with the specified family-wise probability of coverage. 
beta_s.HSD <- TukeyHSD(beta_s, ordered = TRUE)
plot(beta_s.HSD)

#### Stil ordination but let's try and subset this by Site to see effects across time within a same site
#### Site-specific ordination ####
ps.b.r.pcoa.site.plot <- plot_ordination(ps.b.r, ps.b.r.pcoa, axes = c(1,2),
                                         color="Site.name", shape="Month") +
  geom_polygon(aes(fill=Site.name), alpha=0.5) +
  geom_point(alpha=0.5, size=4) +
  facet_wrap(~Site.name, nrow = 3, ncol =5) +
  scale_alpha(guide = 'none') +
  theme(legend.position = "none")
ps.b.r.pcoa.site.plot
ggsave("SGG-PCoA-Axes12-R-Month_variation_Site_Facets.plot.png", 
       height=16, 
       width=20)

#### Unconstrained ordination ####
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
#RDA
sgg_metadata_sig<-read.csv("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_metadata/SGG_metadata_sig.csv", header=TRUE)
#reorder sgg_metadata based on seqtab_bimR
order<-row.names(SGG_SVt)
sgg_metadata_sig_ord<-sgg_metadata_sig[match(order, sgg_metadata_sig$Sample_Code),]
rownames(sgg_metadata_sig_ord) <- sgg_metadata_sig_ord$Sample_Code

ps_sig <- phyloseq(otu_table(SGG_SVt, taxa_are_rows=FALSE),
                   sample_data(sgg_metadata_sig_ord), 
                   tax_table(SGG_Tax))
ps_sig

sample_data(ps_sig)$Month = factor(sample_data(ps_sig)$Month, 
                                   levels = c("April","August","December","Control"))
ps_sig_ntf = subset_samples(ps_sig, Site.name != "Taffs Well PUMPED")
ps_sig_ntf_nc = subset_samples(ps_sig_ntf, Sample_type != "Control")
ps_sig_ntf_nc.r = transform_sample_counts(ps_sig_ntf_nc, function(x) x/sum(x))

#CAP
CAP_sig_ps.b.r <- ordinate(ps_sig_ntf_nc.r,
                           method = "CAP",
                           distance = "jsd",
                           k=3, trymax=1e3, weighted=TRUE,
                           formula = ~ Site.name + Month + Temperature + DO + pH + EC + Cl + SO4 + X6.Li)

CAP_sig_plot_a12 <- plot_ordination(ps_sig_ntf_nc.r,
                                    ordination = CAP_sig_ps.b.r, color = "Site.name", axes = c(1,2)) +
  aes(shape = Month) +
  geom_point(aes(colour = Site.name), alpha = 0.4, size = 4) +
  geom_polygon(aes(fill=Site.name), alpha=0.5) +
  geom_point(aes(colour = Site.name), size = 1.5)

arrowmat_CAP <- vegan::scores(CAP_sig_ps.b.r, display = "bp")
arrowmat_CAP<-arrowmat_CAP[!grepl("Site.name", rownames(arrowmat_CAP)),]
arrowmat_CAP<-arrowmat_CAP[!grepl("Month", rownames(arrowmat_CAP)),]
arrowdf_CAP <- data.frame(labels = rownames(arrowmat_CAP), arrowmat_CAP)
arrow_map_CAP <- aes(xend = CAP1, yend = CAP2,
                     x = 0, y = 0,
                     shape = NULL, color = NULL,
                     label = labels)
label_map_CAP <- aes(x = 1.3 * CAP1, y = 1.3 * CAP2,
                     shape = NULL, color = NULL,
                     label = labels)
arrowhead_CAP = arrow(length = unit(0.02, "npc"))

CAP_sig_plot_a12 = CAP_sig_plot_a12 +
  geom_segment(mapping = arrow_map_CAP, size = .5, 
               data = arrowdf_CAP, color = "gray", arrow = arrowhead_CAP) +
  geom_text(mapping = label_map_CAP, size = 4, 
            data = arrowdf_CAP, show.legend = FALSE)

arrowmat_CAP_13 <- vegan::scores(CAP_sig_ps.b.r, display = "bp", choices=c(1,3))
arrowmat_CAP_13<-arrowmat_CAP_13[!grepl("Site.name", rownames(arrowmat_CAP_13)),]
arrowmat_CAP_13<-arrowmat_CAP_13[!grepl("Month", rownames(arrowmat_CAP_13)),]
arrowdf_CAP_13 <- data.frame(labels = rownames(arrowmat_CAP_13), arrowmat_CAP_13)
arrow_map_CAP_13 <- aes(xend = CAP1, yend = CAP3,
                        x = 0, y = 0,
                        shape = NULL, color = NULL,
                        label = labels)
label_map_CAP_13 <- aes(x = 1.3 * CAP1, y = 1.3 * CAP3,
                        shape = NULL, color = NULL,
                        label = labels)
arrowhead_CAP_13 = arrow(length = unit(0.02, "npc"))
CAP_sig_plot_a13 <- plot_ordination(ps_sig_ntf_nc.r,  ordination = CAP_sig_ps.b.r, color = "Site.name",axes = c(1,3)) +
  aes(shape = Month) +
  geom_point(aes(colour = Site.name), alpha = 0.4, size = 4) +
  geom_polygon(aes(fill=Site.name), alpha=0.5) +
  geom_point(aes(colour = Site.name), size = 1.5) +
  geom_segment(mapping = arrow_map_CAP_13, size = .5,
               data = arrowdf_CAP_13, color = "gray", arrow = arrowhead_CAP_13) +
  geom_text(mapping = label_map_CAP_13, size = 4,
            data = arrowdf_CAP_13, show.legend = FALSE)
ps.b.r.CAP_sig.plot<-plot_grid(CAP_sig_plot_a12, CAP_sig_plot_a13,
                               labels=c("A","B"))
ps.b.r.CAP_sig.plot
ggsave("SGG-CAP-EnvVars-Axes123-R-Month+Site.plot.png",
       height=8,
       width=20)

# 4 - Temperature effects on ps? ----
ps.b.temp<-ps.b
sample_data(ps.b.temp)$Site.name = factor(sample_data(ps.b.temp)$Site.name, 
                                          levels = c("Blaenavon","Cefn Hengoed","Mountain Gate","Taff Bargoed",
                                                     "Glyncastle","Dinas","Ynysarwed","Celynen North","Morlais","Lindsay",
                                                     "Six Bells","Crumlin Navigation","Taffs Well"))
temp_table<-sample_data(ps.b.temp)[,c("Sample_Code","Month","Site.name","Temp.C")]
temp_table.plot<-ggplot(temp_table, 
                        aes(x=Site.name, y=Temp.C, fill=Month)) +
  geom_point(aes(fill=Month, color=Month),
             position = position_dodge(width=0.5)) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1)) +
  geom_hline(yintercept = 13.4,
             color="dark red",
             linetype="dashed") +
  geom_hline(yintercept = 11,
             color="gray",
             linetype="dashed") +
  annotate("text", "Mountain Gate", 13.5, vjust = -1, hjust=.45, 
           label = "Average mine water temperature (Farr et al, 2016)") +
  labs(y=expression("Temperature in"*~degree*C)) +
  ggtitle("Temperatures across the SGG sampled sites")
temp_table.plot
ggsave("Temperatures across the SGG sampled sites.png")

#Relating SVs to temp ====
#Adonis to relate SVs to env variables
sgg_metadata_adonis<-read.csv("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_metadata/SGG_metadata_adonis.csv", header=TRUE)
#reorder sgg_metadata based on seqtab_bimR
order<-row.names(SGG_SVt)
sgg_metadata_adonis_ord<-sgg_metadata_adonis[match(order, sgg_metadata_adonis$Sample_Code),]
rownames(sgg_metadata_adonis_ord) <- sgg_metadata_adonis_ord$Sample_Code

ps_adonis <- phyloseq(otu_table(SGG_SVt, taxa_are_rows=FALSE),
                      sample_data(sgg_metadata_adonis_ord), 
                      tax_table(SGG_Tax))
ps_adonis

sample_data(ps_adonis)$Month = factor(sample_data(ps_adonis)$Month, 
                                      levels = c("April","August","December","Control"))
ps_adonis_ntf = subset_samples(ps_adonis, Site.name != "Taffs Well PUMPED")
ps_adonis_ntf_nc = subset_samples(ps_adonis_ntf, Sample_type != "Control")

ps_adonis_ntf_nc.r = transform_sample_counts(ps_adonis_ntf_nc, function(x) x/sum(x))
df = as(sample_data(ps_adonis_ntf_nc.r), "data.frame")
df$Temp.C<-as.numeric(df$Temp.C)
d = phyloseq::distance(ps_adonis_ntf_nc, "bray")
ps.adonis.ntf.nc.adonis = adonis(d ~ Month + Site.name + Temp.C + pH + EC.uS.cm + DO..sat +K.µg.ml +Na.µg.ml+Cl..µg.ml+SO4..µg.ml+X6.Li.ng.ml+X7.Li.ng.ml+X7.Li..He..ng.ml+X11.B.ng.ml, 
                                 df, strata=df$Time_point)

ps.adonis.ntf.nc.adonis
# plot(ps.adonis.ntf.nc.adonis$aov.tab)
# densityplot(permustats(ps.adonis.ntf.nc.adonis))

#Relating SVs to temp p.2 ====
#now from https://rpubs.com/dillmcfarlan/R_microbiotaSOP
#see what SVs correlate to temperature with cor()
SV1 = as(otu_table(ps_adonis_ntf_nc.r), "matrix")
# transpose if necessary
if(taxa_are_rows(ps_adonis_ntf_nc.r)){SV1 <- t(SV1)}
# Coerce to data.frame
SV1df = as.data.frame(SV1)
#Kendall's Tau doesn't assume normality - ideal for microbioshit
cor.kendall.ps_adonis_ntf_nc.r = cor(SV1df, 
                                     sample_data(ps_adonis_ntf_nc.r)$Temp.C, 
                                     method = "kendall")
# plot(sort(cor.kendall.ps_adonis_ntf_nc.r, decreasing = TRUE))
# abline(h = c(-0.5,0, 0.5))
cor.kendall.ps_adonis_ntf_nc.r.df<-as.data.frame(cor.kendall.ps_adonis_ntf_nc.r)
cor.kendall.ps_adonis_ntf_nc.r.df.ord<-cor.kendall.ps_adonis_ntf_nc.r.df[order(cor.kendall.ps_adonis_ntf_nc.r.df$V1), , drop = FALSE]
cor.kendall.ps_adonis_ntf_nc.r.df.ord.InterestSVs<-head(cor.kendall.ps_adonis_ntf_nc.r.df.ord, n=5)

kendall.SVs<-rownames(cor.kendall.ps_adonis_ntf_nc.r.df.ord.InterestSVs)
kendall.SVs.ps.b.r = prune_taxa(kendall.SVs, ps.b.r)
tax_table(kendall.SVs.ps.b.r)
#The 5 Svs are: n=3 Gallionella, n=1 Nitrospira, n=1 Omnitrophica (NA at Genus level)
cors <- ddply(mtcars, c("vs", "am"), summarise, cor = round(cor(mpg, wt), 2))

#just for plotting purposes, change name of SVs
taxa_names(kendall.SVs.ps.b.r) <- paste0("SV_", seq(ntaxa(kendall.SVs.ps.b.r)))
# kendall.SVs.ps.b.r.glom <- tax_glom(kendall.SVs.ps.b.r, taxrank = 'Species', NArm = FALSE)
kendall.SVs.ps.b.r.glom.psdf <- data.table(psmelt(kendall.SVs.ps.b.r))
kendall.SVs.ps.b.r.glom.psdf$Species <- as.character(kendall.SVs.ps.b.r.glom.psdf$Species)

kendall.SVs.ps.b.r.glom.psdf.plot<-ggplot(kendall.SVs.ps.b.r.glom.psdf[Abundance > 0], 
                                          aes(x = Temp.C, y = Abundance, fill=OTU)) + 
  # facet_grid(~Month, scales="free_x") +
  geom_point(aes(color=OTU, fill=OTU,
                 size=Abundance)) +
  geom_smooth(aes(color=OTU, fill=OTU),
              method='lm', se=FALSE) +
  theme(axis.title.x = element_blank())+ 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  scale_y_log10() +
  ylab("Relative Abundance\n")

#create labels with Kendall's correlation
#on hold, waiting for reply from phyloseq people

kendall.SVs.ps.b.r.glom.psdf.plot
# annotate("text", x=14, y=0.02, label = "R^2=0.78") +
# annotate("text", x=14, y=0.01, label = "alpha=0.00") +
# annotate("text", x=14, y=0.03, label = "beta=0.67")

ggsave("SGG-KendallsCorSVs-TempOnly.plot.png", 
       height=8, 
       width=20)

#Relating SVs to Cl ====
#now from https://rpubs.com/dillmcfarlan/R_microbiotaSOP
cor.kendall.Cl.ps_adonis_ntf_nc.r = cor(SV1df, 
                                        sample_data(ps_adonis_ntf_nc.r)$Cl..µg.ml, 
                                        method = "kendall")
plot(sort(cor.kendall.Cl.ps_adonis_ntf_nc.r, decreasing = TRUE))
abline(h = c(-0.5,0, 0.5))
cor.kendall.Cl.ps_adonis_ntf_nc.r.df<-as.data.frame(cor.kendall.Cl.ps_adonis_ntf_nc.r)
cor.kendall.Cl.ps_adonis_ntf_nc.r.df.ord<-cor.kendall.Cl.ps_adonis_ntf_nc.r.df[order(cor.kendall.Cl.ps_adonis_ntf_nc.r.df$V1), , 
                                                                               drop = FALSE]
cor.kendall.Cl.ps_adonis_ntf_nc.r.df.ord.InterestSVs<-head(cor.kendall.Cl.ps_adonis_ntf_nc.r.df.ord, n=22)
#22 SVs >-0.5
kendall.Cl.SVs<-rownames(cor.kendall.Cl.ps_adonis_ntf_nc.r.df.ord.InterestSVs)
kendall.Cl.SVs.ps.b.r = prune_taxa(kendall.Cl.SVs, ps.b.r)
as.data.frame(tax_table(kendall.Cl.SVs.ps.b.r))$Phylum
#The 22 Svs are: 1 Proteobacteria, 10 Omnitrophica, 6 Parcubacteria, 2 Nitrospirae, 1 Gracillibacteria, 1 WS2, 1 Ignavibacteriae

#just for plotting purposes, change name of SVs
taxa_names(kendall.Cl.SVs.ps.b.r) <- paste0("SV_", seq(ntaxa(kendall.Cl.SVs.ps.b.r)))
# kendall.SVs.ps.b.r.glom <- tax_glom(kendall.SVs.ps.b.r, taxrank = 'Species', NArm = FALSE)
kendall.Cl.SVs.ps.b.r.glom.psdf <- data.table(psmelt(kendall.Cl.SVs.ps.b.r))
kendall.Cl.SVs.ps.b.r.glom.psdf$Species <- as.character(kendall.Cl.SVs.ps.b.r.glom.psdf$Species)
kendall.Cl.SVs.ps.b.r.glom.psdf.plot<-ggplot(kendall.Cl.SVs.ps.b.r.glom.psdf[Abundance > 0], 
                                             aes(x = Cl..µg.ml, y = Abundance, fill=OTU)) + 
  geom_point(aes(color=OTU, fill=OTU,
                 size=Abundance)) +
  geom_smooth(aes(color=OTU, fill=OTU),
              method='lm', se=FALSE) +
  theme(axis.title.x = element_blank())+ 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  scale_y_log10() +
  ylab("Relative Abundance\n")

# kendall.Cl.SVs.ps.b.r.glom.psdf.plot

# Prevalence in the dataset ----
#Completely copied from http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html#filter-taxa
mdt = fast_melt(ps.b)
prevdt = mdt[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
#Test
# ggplot(prevdt, aes(Prevalence)) + 
#   geom_histogram() + 
#   ggtitle("Histogram of Taxa Prevalence")
prevdt[(Prevalence <= 0), .N]
prevdt[(Prevalence <= 1), .N]
prevdt[(Prevalence <= 2), .N]

# taxa cumulative sum
prevcumsum = prevdt[, .N, by = Prevalence]
setkey(prevcumsum, Prevalence)
prevcumsum[, CumSum := cumsum(N)]
pPrevCumSum = ggplot(prevcumsum, aes(Prevalence, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Prevalence") +
  ylab("SVs Filtered") +
  ggtitle("SVs that would be filtered vs. the minimum count threshold")
pPrevCumSum
#
addPhylum = unique(copy(mdt[, list(TaxaID, Phylum)]))
# Join by TaxaID
setkey(prevdt, TaxaID)
setkey(addPhylum, TaxaID)
prevdt <- addPhylum[prevdt]
showPhyla = prevdt[, sum(TotalCounts), by = Phylum][order(-V1)][1:9]$Phylum
setkey(prevdt, Phylum)
prev_totcounts<-ggplot(prevdt[showPhyla], 
                       mapping = aes(Prevalence, TotalCounts, color = Phylum)) + 
  geom_point(size = 4, alpha = 0.7) + 
  scale_y_log10()
prev_totcounts
ggsave("SGG-PrevalenceVSTotalCounts.plot.png", 
       height=8, 
       width=20)

#### Networking! ####
ps.b.r.f = filter_taxa(ps.b.r, function(x) sum(x) > .01, TRUE)
taxa_names(ps.b.r.f) <- paste("SV-", 1:ntaxa(ps.b.r.f), sep="")
ps.b.r.f.jsd <- phyloseq::distance(ps.b.r.f, method = "jaccard")


#save this for when issue https://github.com/joey711/phyloseq/issues/885 is answered
taxa_net<-plot_net(ps.b.r.f, maxdist = 0.8, type="taxa", 
         laymeth="fruchterman.reingold",
         rescale = TRUE,
         point_alpha = 0.6)
taxa_net

#samples

sample_net.ig <- make_network(ps.b.r.f, max.dist=0.8)

V(sample_net.ig)$Site.name=as.character(sample_data(ps.b.r.f)$Site.name[match(V(sample_net.ig)$name,
                                  sample_data(ps.b.r.f)$Sample_Code)])
V(sample_net.ig)$Month=as.character(sample_data(ps.b.r.f)$Month[match(V(sample_net.ig)$name,
                                  sample_data(ps.b.r.f)$Sample_Code)])
#this is probably wrong, but will do for now (?)
E(sample_net.ig)$weight <- strength(sample_net.ig, mode = "all")
E(sample_net.ig)$weight <- E(sample_net.ig)$weight/18

sample_net.ig.plot<-ggnet2(sample_net.ig, 
                           color="Site.name", shape="Month", 
                           size=5, alpha = 0.8, 
                           edge.size = "weight") +
    scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(13)) +
  theme(legend.text = element_text(size=10))
sample_net.ig.plot

#### Temporal trends ####
#first indicator of temporal trends might be alpha-diversity
plot_richness(ps.b, x="Month", measures="Observed",
              color="Site.name") + 
  facet_wrap(~Site.name) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text("Shannon Diversity"))

plot_richness(ps.b, x="Month", measures="Shannon",
              color="Site.name") + 
  geom_point(size = 5, alpha=0.4) +
  facet_wrap(~Site.name) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1),
        axis.title.x = element_blank()) +
  ylab("Shannon Diversity")
ggsave("SGG-ShannonDiv-Month.plot.png", 
       height=8, 
       width=16)

plot_richness(ps.b, x="Month", measures="Simpson",
              color="Site.name") + 
  facet_wrap(~Site.name) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, vjust=1, hjust=1),
        axis.title.x = element_blank(),
        axis.title.y = element_text("Shannon Diversity"))

#It seems that alpha diversity does change across time in some sites
#how about beta diversity?

# Join sample data and ordination axes together in one data.table
ps.b.r.pcoa.v = data.table(ps.b.r.pcoa$vectors, keep.rownames = TRUE)
sdt = data.table(as(sample_data(ps.b.r), "data.frame"), keep.rownames = TRUE)

setnames(ps.b.r.pcoa.v, "rn", "Sample.ID")
df <- subset(sdt, select = -c(Sample.ID) )
setnames(df, "rn", "Sample.ID")
setkey(ps.b.r.pcoa.v, Sample.ID)
setkey(df, Sample.ID)
ps.b.r.pcoa.v <- ps.b.r.pcoa.v[df]
# setorder(ps.b.r.pcoa.v, Month)
# Axis 1
ggplot(ps.b.r.pcoa.v, aes(Month, Axis.1, color = Site.name)) +
  geom_point(size = 3) +
  geom_point(size = 5, alpha=0.4) + 
  facet_wrap(~Site.name) +
  theme(axis.text.x = element_text(angle = 25, vjust=1, hjust=1)) +
  ylab("JSD PCoA Axis 1")
ggsave("SGG-JSDPCoAAxis1-Month.plot.png", 
       height=8, 
       width=16)

#can we detect differential SV abundance/clade relative abundance across time with DeSeq2?
month_ps_dds = phyloseq_to_deseq2(ps.b, ~ Month)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(month_ps_dds), 1, gm_mean)
month_ps_dds = estimateSizeFactors(month_ps_dds, geoMeans = geoMeans)
month_ps_dds = DESeq(month_ps_dds, fitType="local")

res = results(month_ps_dds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps.b)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Order order
x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))
ggplot(sigtab, aes(x=Order, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

#only 4 SVs significantly differ across month of sampling.
#quick track down of these nasty buggers
month_sig_SVs <- rownames(sigtab)
ex2 <- prune_taxa(month_sig_SVs, ps.b.r)

plot_bar(ex2, x="Site.name", fill="Phylum") +
  facet_wrap(~Month) +
  theme(axis.text.x = element_text(angle = 35, vjust=1, hjust=1),
        axis.title.x = element_blank())
ggsave("SGG-DiffAbDESeq2SVs-Month.plot.png", 
       height=8, 
       width=16)
