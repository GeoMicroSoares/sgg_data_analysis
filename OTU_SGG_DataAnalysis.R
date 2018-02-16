#completely adapted from https_otu://benjjneb.github.io/dada2/tutorial.html
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

setwd("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_OTU")

SGG_otu_tab = "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_OTU/SGG_OTU_m.txt"
SGG_otu_tab_imp<-read.cotu(SGG_otu_tab, sep="\t", header = TRUE, row.names = 1, check.names = FALSE)
OTU = otu_table(SGG_otu_tab_imp, taxa_are_rows = TRUE)
OTU_t = t(OTU)
SGG_otu_tax = "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_OTU/SGG_OTU_tax_c.txt"
SGG_otu_tax_imp<-read.cotu(SGG_otu_tax, sep="\t", header = TRUE, row.names = 1)
SGG_OTU_tax<-as.matrix(SGG_otu_tax_imp)

sgg_metadata<-read.cotu("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_metadata/SGG_metadata.cotu", header=TRUE)
order<-rownames(OTU_t)
sgg_metadata_ord<-sgg_metadata[match(order, sgg_metadata$Sample_Code),]
rownames(sgg_metadata_ord) <- sgg_metadata_ord$Sample_Code

ps_otu<-phyloseq(otu_table(OTU), 
                 tax_table(SGG_OTU_tax), 
                 sample_data(sgg_metadata_ord))
ps_otu

# 0 - Basic stats ----
#how many reads left in seqtab_bimR?
sum(rowSums(otu_table(ps_otu)))
#How many OTUs in ps_otu?
ntaxa(otu_table(ps_otu))
#70,937 SVs to 9,656,492 reads

#Seq Depth
ps_otudt = data.table(as(sample_data(ps_otu), "data.frame"),
                  TotalReads = sample_sums(ps_otu), keep.rownames = TRUE)
setnames(ps_otudt, "rn", "SampleID")
ps_otueqDepth = ggplot(ps_otudt, aes(TotalReads)) + 
  geom_histogram(binwidth = 4000) + 
  ggtitle("Sequencing Depth")
ps_otueqDepth
ggsave("ps_otueqDepth.png", width=14, height=6)

#Number of SVs per phylum
sv.phyl<-as.data.frame(table(tax_table(ps_otu)[, "Phylum"], exclude = NULL))
sv.phyl.ord <- sv.phyl[order(-sv.phyl$Freq),] 
sv.phyl.ord

#Number of sequences per SV, Sequences per sample
readsumsdf_otu = data.frame(nreads = sort(taxa_sums(ps_otu), TRUE), sorted = 1:ntaxa(ps_otu), type = "OTUs")
readsumsdf_otu = rbind(readsumsdf_otu, data.frame(nreads = sort(sample_sums(ps_otu), 
                                                        TRUE), sorted = 1:nsamples(ps_otu), type = "Samples"))
p_read_sums = ggplot(readsumsdf_otu, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p_read_sums + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
ggsave("SVsSeq_SeqsSpl.png", width=14, height=6)

#0.5 - Custom datasets  ----
#now without Controls and Taffs Well PUMPED
ps_otu.a = subset_samples(ps_otu, Site.name != "Taffs Well PUMPED")
ps_otu.b = subset_samples(ps_otu.a, Sample_type != "Control")
ps_otu.b <- prune_taxa(taxa_sums(ps_otu.b) > 0, ps_otu.b)
#order months
sample_data(ps_otu.b)$Month = factor(sample_data(ps_otu.b)$Month, 
                                 levels = c("April","August","December"))
#to relative abundances
ps_otu.b.r <-  transform_sample_counts(ps_otu.b, function(x) {x/sum(x)} ) 

# 1 - Taxonomic analysis ----
# Phylum-level tax, transformed ====
# agglomerate taxa
ps_otu.b.r.glom <- tax_glom(ps_otu.b.r, taxrank = 'Phylum', NArm=FALSE)
# create dataframe from phyloseq object
ps_otu.b.r.glom.ps_otudf_otu <- data.table(psmelt(ps_otu.b.r.glom))
# convert Phylum to a character vector from a factor because R
ps_otu.b.r.glom.ps_otudf_otu$Phylum <- as.character(ps_otu.b.r.glom.ps_otudf_otu$Phylum)
# group dataframe by Phylum, calculate median rel. abundance
ps_otu.b.r.glom.ps_otudf_otu[, median := median(Abundance, na.rm = TRUE), 
                 by = "Phylum"]
# Change name to remainder of Phylum less than 1%
ps_otu.b.r.glom.ps_otudf_otu[(median <= 0.01), Phylum := "Other"]

ps_otu.b.r.glom.ps_otudf_otu.plot<-ggplot(ps_otu.b.r.glom.ps_otudf_otu[Abundance > 0], 
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
ps_otu.b.r.glom.ps_otudf_otu.plot
ggsave("sgg_otu-0.01pc-Phylum-R-Month+Site.plot.png", height=8, width=16)

# Phylum-, Proteobacterial class-level tax, not transformed ====
#Removing Proteobacteria from previous table
ps_otu.b.r.glom.ps_otudf_otu.noProteo<-ps_otu.b.r.glom.ps_otudf_otu[!grepl("Proteobacteria", ps_otu.b.r.glom.ps_otudf_otu$Phylum),]
#New table, with Proteobacterial classes only
ps_otu.b.r.ProtCl = subset_taxa(ps_otu.b.r, Phylum == "Proteobacteria")
ps_otu.b.r.ProtCl.glom <- tax_glom(ps_otu.b.r.ProtCl, taxrank = 'Class', NArm=FALSE)
# create dataframe from phyloseq object
ps_otu.b.r.ProtCl.glom.ps_otudf_otu <- data.table(psmelt(ps_otu.b.r.ProtCl.glom))
# convert Class to a character vector from a factor because R
ps_otu.b.r.ProtCl.glom.ps_otudf_otu$Class <- as.character(ps_otu.b.r.ProtCl.glom.ps_otudf_otu$Class)
# group dataframe by Phylum, calculate median rel. abundance
ps_otu.b.r.ProtCl.glom.ps_otudf_otu[, median := median(Abundance, na.rm = TRUE), 
                        by = "Class"]
# Change name to remainder of Phylum less than 1%
ps_otu.b.r.ProtCl.glom.ps_otudf_otu[(median <= 0.01), Class := "Other Proteobacteria"]
#Check if they're all there
# ggplot(ps_otu.b.r.ProtCl.glom.ps_otudf_otu[Abundance > 0], 
#        aes(x = Site.name, y = Abundance/3, fill = Class)) +
#   facet_grid(Month~.) +
#   geom_bar(stat = "identity")
#join the two data frames
ps_otu.b.r.ProtCl.glom.ps_otudf_otu.noP <- subset(ps_otu.b.r.ProtCl.glom.ps_otudf_otu, select = -Phylum)
#same name= allows merging
names(ps_otu.b.r.glom.ps_otudf_otu.noProteo)[names(ps_otu.b.r.glom.ps_otudf_otu.noProteo) == 'Phylum'] <- 'Taxa'
names(ps_otu.b.r.ProtCl.glom.ps_otudf_otu.noP)[names(ps_otu.b.r.ProtCl.glom.ps_otudf_otu.noP) == 'Class'] <- 'Taxa'
#join the two
ps_otu.b.r.ProteoCl.AllPhy <- rbind(ps_otu.b.r.glom.ps_otudf_otu.noProteo, ps_otu.b.r.ProtCl.glom.ps_otudf_otu.noP)

ps_otu.b.r.ProteoCl.AllPhy.plot<-ggplot(ps_otu.b.r.ProteoCl.AllPhy[Abundance > 0], 
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
ps_otu.b.r.ProteoCl.AllPhy.plot
ggsave("sgg_otu-0.01pc-Phylum+ProtoCl-R-Month+Site.plot.png", height=8, width=16)

# Genus-level tax, transformed ====
#Tracking taxa responsible for Beta, Eps_otuilonproteobacteria dominances
# Gallionellaceae, not transformed ####
ps_otu.b.gallionelac = subset_taxa(ps_otu.b, Family=="Gallionellaceae")
ps_otu.b.gallionelac.plot<-plot_bar(ps_otu.b.gallionelac, 
                                fill="Genus", facet_grid=~Month, x="Site.name") + 
  geom_bar(aes(color=Genus, fill=Genus), 
           stat="identity", position="stack")
ps_otu.b.gallionelac.plot

# Campylobacterales, not transformed ####
ps_otu.b.helicob = subset_taxa(ps_otu.b, Family=="Helicobacteraceae")
ps_otu.b.helicob.plot<-plot_bar(ps_otu.b.helicob, 
                            fill="Genus", facet_grid=~Month, x="Site.name") + 
  geom_bar(aes(color=Genus, fill=Genus), 
           stat="identity", position="stack")
ps_otu.b.helicob.plot

ps_otu.b.helicobsv.phyl<-as.data.frame(table(tax_table(ps_otu.b.helicob)[, "Genus"], exclude = NULL))
ps_otu.b.helicobsv.phyl.ord <- ps_otu.b.helicobsv.phyl[order(-ps_otu.b.helicobsv.phyl$Freq),]
ps_otu.b.helicobsv.phyl.ord

#How about all these together?
ps_otu.b.r.HelicobGall = subset_taxa(ps_otu.b.r, Family=="Helicobacteraceae" | Family=="Gallionellaceae")
ps_otu.b.r.HelicobGall.glom <- tax_glom(ps_otu.b.r.HelicobGall, taxrank = 'Genus', NArm=FALSE)
ps_otu.b.r.HelicobGall.glom.ps_otudf_otu <- data.table(psmelt(ps_otu.b.r.HelicobGall.glom))
ps_otu.b.r.HelicobGall.glom.ps_otudf_otu$Genus <- as.character(ps_otu.b.r.HelicobGall.glom.ps_otudf_otu$Genus)

ps_otu.b.r.HelicobGall.glom.ps_otudf_otu.plot<-ggplot(ps_otu.b.r.HelicobGall.glom.ps_otudf_otu,
                                          aes(x=Genus,y=Abundance,fill=Genus)) + 
  geom_point(aes(color=Genus, fill=Genus,
                 size=Abundance)) + 
  facet_grid(Month~Site.name)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ps_otu.b.r.HelicobGall.glom.ps_otudf_otu.plot
ggsave("sgg_otu-Helicobacteraceae+Gallionellaceae-Genera-R-Month-Site.plot.png",
       width = 18,
       height = 8)

# Proteobacteria-less tax, transformed ====
ps_otu.b.r.noProt<-subset_taxa(ps_otu.b.r, Phylum!="Proteobacteria")
#check phyla
sv.phyl.noProt<-as.data.frame(table(tax_table(ps_otu.b.r.noProt)[, "Phylum"], exclude = NULL))
sv.phyl.noProt.ord <- sv.phyl.noProt[order(-sv.phyl.noProt$Freq),] 
#let's go for Phyla w/more than 1000 SV's following Proteobacteria
ps_otu.b.r.noProt.glom <- tax_glom(ps_otu.b.r.noProt, taxrank = 'Phylum')
ps_otu.b.r.noProt.glom.ps_otudf_otu <- data.table(psmelt(ps_otu.b.r.noProt.glom))
ps_otu.b.r.noProt.glom.ps_otudf_otu$Phylum <- as.character(ps_otu.b.r.noProt.glom.ps_otudf_otu$Phylum)
ps_otu.b.r.noProt.glom.ps_otudf_otu[, median := median(Abundance, na.rm = TRUE), 
                        by = "Phylum"]
ps_otu.b.r.noProt.glom.ps_otudf_otu[(median <= 0.005), Phylum := "Other"]

ps_otu.b.r.noProt.glom.ps_otudf_otu.plot<-ggplot(ps_otu.b.r.noProt.glom.ps_otudf_otu[Abundance > 0], 
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
ps_otu.b.r.noProt.glom.ps_otudf_otu.plot
ggsave("sgg_otu-PhylaBelowProteobacteria-R-Month+Site.plot.png", 
       height=8, 
       width=20)
#
# Archaea, non-transformed ====
ps_otu.b.r.archaea <- subset_taxa(ps_otu.b.r, Domain=="Archaea")

ps_otu.b.r.archaea.phylum<-as.data.frame(table(tax_table(ps_otu.b.r.archaea)[, "Phylum"], exclude = NULL))
ps_otu.b.r.archaea.phylum.ord <- ps_otu.b.r.archaea.phylum[order(-ps_otu.b.r.archaea.phylum$Freq),] 
ps_otu.b.r.archaea.phylum.ord

ps_otu.b.r.archaea.glom <- tax_glom(ps_otu.b.r.archaea, taxrank = 'Phylum', NArm = FALSE )
ps_otu.b.r.archaea.glom.ps_otudf_otu <- data.table(psmelt(ps_otu.b.r.archaea.glom))
ps_otu.b.r.archaea.glom.ps_otudf_otu$Phylum <- as.character(ps_otu.b.r.archaea.glom.ps_otudf_otu$Phylum)
# plot(sort(ps_otu.b.r.OmniParcu.glom.ps_otudf_otu$Abundance, decreasing = TRUE))
#plot only classes representing >1% of the sample profile in the future
ps_otu.b.r.archaea.glom.ps_otudf_otu.plot<-ggplot(ps_otu.b.r.archaea.glom.ps_otudf_otu[Abundance > 0], 
                                      aes(x = Phylum, y = Abundance, fill = Phylum)) + 
  facet_grid(Month~Site.name, space = "free_x") +
  geom_point(aes(color=Phylum, fill=Phylum,
                 size=Abundance)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+ 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance\n")+
  geom_hline(yintercept = 0.01,
             color="black",
             linetype="dashed")
#totals vary because the bars = medians
ps_otu.b.r.archaea.glom.ps_otudf_otu.plot
ggsave("sgg_otu-ArchaeaClasses-R-Month+Site.plot.png", 
       height=8, 
       width=20)

# Nitrospirae, non-transformed ====
ps_otu.b.r.nitrospirae <- subset_taxa(ps_otu.b.r, Phylum=="Nitrospirae")
#
ps_otu.b.r.nitrospirae.Family<-as.data.frame(table(tax_table(ps_otu.b.r.nitrospirae)[, "Family"], exclude = NULL))
ps_otu.b.r.nitrospirae.Family.ord <- ps_otu.b.r.nitrospirae.Family[order(-ps_otu.b.r.nitrospirae.Family$Freq),] 
ps_otu.b.r.nitrospirae.Family.ord
#
ps_otu.b.r.nitrospirae.glom <- tax_glom(ps_otu.b.r.nitrospirae, taxrank = 'Family', NArm = FALSE )
ps_otu.b.r.nitrospirae.glom.ps_otudf_otu <- data.table(psmelt(ps_otu.b.r.nitrospirae.glom))
ps_otu.b.r.nitrospirae.glom.ps_otudf_otu$Family <- as.character(ps_otu.b.r.nitrospirae.glom.ps_otudf_otu$Family)
#plot only Familyes representing >1% of the sample profile in the future
ps_otu.b.r.nitrospirae.glom.ps_otudf_otu.plot<-ggplot(ps_otu.b.r.nitrospirae.glom.ps_otudf_otu[Abundance > 0], 
                                          aes(x = Family, y = Abundance, fill = Family)) + 
  facet_grid(Month~Site.name, space = "free_x") +
  geom_point(aes(color=Family, fill=Family,
                 size=Abundance)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+ 
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance\n")+
  geom_hline(yintercept = 0.01,
             color="black",
             linetype="dashed")
#totals vary because the bars = medians
ps_otu.b.r.nitrospirae.glom.ps_otudf_otu.plot
ggsave("sgg_otu-NitrospiraeFamily-R-Month+Site.plot.png", 
       height=8, 
       width=20)

# 3 - Ordination on ps_otu ----
# nMDS and PCoA on dataset
# ps_otu.b.r.dist <- phyloseq::distance(ps_otu.b.r, method="bray")
# Calculate ordination
ps_otu.b.r.pcoa  <- ordinate(ps_otu.b.r, "PCoA", distance="jsd", 
                         k=3, trymax=1e3, weighted=TRUE)
#How good is this, i.e. how much of the variance is explained by the first few components?
plot_scree(ps_otu.b.r.pcoa)
evals<-ps_otu.b.r.pcoa$values$Eigenvalues
#Axis 1 and 2
pcoa.a12<-sum(evals[1:2])/sum(evals)
pcoa.a12.r<-round(pcoa.a12, digits=4)*100
#Axis 1 and 3
pcoa.a13<-sum(evals[c(1,3)])/sum(evals)
pcoa.a13.r<-round(pcoa.a13, digits=4)*100

ps_otu.b.r.pcoa.plot.a1 <- plot_ordination(ps_otu.b.r, ps_otu.b.r.pcoa, axes = c(1,2),
                                       color="Site.name", shape="Month") +
  geom_polygon(aes(fill=Site.name), alpha=0.5) +
  geom_point(alpha=0.5, size=4)+
  scale_alpha(guide = 'none') +
  annotate("text", y= 0.3, x =-0.3,label=paste("Axes 1 and 2:\n",pcoa.a12.r,"% variation")
           ,hjust=0.5, size=5)

ps_otu.b.r.pcoa.plot.a2 <- plot_ordination(ps_otu.b.r, ps_otu.b.r.pcoa, axes = c(1,3),
                                       color="Site.name", shape="Month")+
  geom_polygon(aes(fill=Site.name), alpha=0.5) +
  geom_point(alpha=0.5, size=4) +
  scale_alpha(guide = 'none')+
  annotate("text", y= 0.15, x =-0.3,label=paste("Axes 1 and 3:\n",pcoa.a13.r,"% variation")
           ,hjust=0.5, size=5)

ps_otu.b.r.pcoa.plot<-plot_grid(ps_otu.b.r.pcoa.plot.a1, ps_otu.b.r.pcoa.plot.a2, 
                            labels=c("A","B"))
ps_otu.b.r.pcoa.plot
ggsave("sgg_otu-PCoA-Axes123-R-Month+Site.plot.png", 
       height=6, 
       width=16)
#
ps_otu.b.r.nmds  <- ordinate(ps_otu.b.r, "NMDS", distance="bray", 
                         k=3, trymax=1e3, sratmax=0.999999)
ps_otu.b.r.nmds.plot.a1 <- plot_ordination(ps_otu.b.r, ps_otu.b.r.nmds, axes = c(1,2),
                                       color="Site.name", shape="Month") +
  geom_polygon(aes(fill=Site.name), alpha=0.5) +
  geom_point(alpha=0.5, size=4)+
  scale_alpha(guide = 'none')
ps_otu.b.r.nmds.plot.a2 <- plot_ordination(ps_otu.b.r, ps_otu.b.r.nmds, axes = c(1,3),
                                       color="Site.name", shape="Month")+
  geom_polygon(aes(fill=Site.name), alpha=0.5) +
  geom_point(alpha=0.5, size=4) +
  scale_alpha(guide = 'none')

ps_otu.b.r.nmds.plot<-plot_grid(ps_otu.b.r.nmds.plot.a1, ps_otu.b.r.nmds.plot.a2, 
                            labels=c("A","B"))
ps_otu.b.r.nmds.plot
ggsave("sgg_otu-nMDS-Axes123-R-Month+Site.plot.png", 
       height=8, 
       width=20)

#see it in 3D, crap foe 
# BCxyz = scores(ps_otu.b.r.nmds, display="sites")
# plot_ly(x=BCxyz[,1], y=BCxyz[,2], z=BCxyz[,3], 
#         type="scatter3d", mode="markers", 
#         color=sample_data(ps_otu.b.r)$Site.name, colors=c("blue","green","red","green",
#                                                       "pink","black","gray","yellow",
#                                                       "orange","lightblue","darkblue","darkgreen",
#                                                       "lightgreen"))

# 3 - Temperature effects on ps_otu? ----
ps_otu.b.temp<-ps_otu.b
sample_data(ps_otu.b.temp)$Site.name = factor(sample_data(ps_otu.b.temp)$Site.name, 
                                          levels = c("Blaenavon","Cefn Hengoed","Mountain Gate","Taff Bargoed",
                                                     "Glyncastle","Dinas","Ynysarwed","Celynen North","Morlais","Lindsay",
                                                     "Six Bells","Crumlin Navigation","Taffs Well"))
temp_table<-sample_data(ps_otu.b.temp)[,c("Sample_Code","Month","Site.name","Temp.C")]
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
  ggtitle("Temperatures across the sgg_otu sampled sites")
temp_table.plot
ggsave("Temperatures across the sgg_otu sampled sites.png")

#Relating SVs to temp ====
#Adonis to relate SVs to env variables
sgg_otu_metadata_adonis<-read.csv("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_metadata/SGG_metadata_adonis.csv",
                                  header=TRUE)
#reorder sgg_otu_metadata based on seqtab_bimR
order<-rownames(OTU_t)
sgg_otu_metadata_adonis_ord<-sgg_otu_metadata_adonis[match(order, sgg_otu_metadata_adonis$Sample_Code),]
rownames(sgg_otu_metadata_adonis_ord) <- sgg_otu_metadata_adonis_ord$Sample_Code

ps_otu_adonis <- phyloseq(otu_table(OTU),
                      sample_data(sgg_otu_metadata_adonis_ord), 
                      tax_table(SGG_OTU_tax))
ps_otu_adonis

sample_data(ps_otu_adonis)$Month = factor(sample_data(ps_otu_adonis)$Month, 
                                      levels = c("April","August","December","Control"))
ps_otu_adonis_ntf = subset_samples(ps_otu_adonis, Site.name != "Taffs Well PUMPED")
ps_otu_adonis_ntf_nc = subset_samples(ps_otu_adonis_ntf, Sample_type != "Control")

ps_otu_adonis_ntf_nc.r = transform_sample_counts(ps_otu_adonis_ntf_nc, function(x) x/sum(x))
df_otu = as(sample_data(ps_otu_adonis_ntf_nc.r), "data.frame")
d = phyloseq::distance(ps_otu_adonis_ntf_nc, "bray")
ps_otu.adonis.ntf.nc.adonis = adonis(d ~ Month + Site.name + Temp.C + pH + EC.uS.cm + DO..sat +K.µg.ml +Na.µg.ml+Cl..µg.ml+SO4..µg.ml+X6.Li.ng.ml+X7.Li.ng.ml+X7.Li..He..ng.ml+X11.B.ng.ml, df_otu)

ps_otu.adonis.ntf.nc.adonis
plot(ps_otu.adonis.ntf.nc.adonis$aov.tab)
densityplot(permustats(ps_otu.adonis.ntf.nc.adonis))

#Relating SVs to temp p.2 ====
#now from https_otu://rpubs.com/dillmcfarlan/R_microbiotaSOP
#see what SVs correlate to temperature with cor()
OTU1 = as(otu_table(ps_otu_adonis_ntf_nc.r), "matrix")
# transpose if necessary
if(taxa_are_rows(ps_otu_adonis_ntf_nc.r)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTU1df_otu = as.data.frame(OTU1)
#Kendall's Tau doesn't assume normality - ideal for microbioshit
cor.kendall.ps_otu_adonis_ntf_nc.r = cor(OTU1df_otu, 
                                     sample_data(ps_otu_adonis_ntf_nc.r)$Temp.C, 
                                     method = "kendall")
plot(sort(cor.kendall.ps_otu_adonis_ntf_nc.r, decreasing = TRUE))
abline(h = c(-0.5,0, 0.5))
cor.kendall.ps_otu_adonis_ntf_nc.r.df_otu<-as.data.frame(cor.kendall.ps_otu_adonis_ntf_nc.r)
cor.kendall.ps_otu_adonis_ntf_nc.r.df_otu.ord<-cor.kendall.ps_otu_adonis_ntf_nc.r.df_otu[order(cor.kendall.ps_otu_adonis_ntf_nc.r.df_otu$V1), , drop = FALSE]
cor.kendall.ps_otu_adonis_ntf_nc.r.df_otu.ord.InterestOTUs<-head(cor.kendall.ps_otu_adonis_ntf_nc.r.df_otu.ord, n=35)

kendall.OTUs<-rownames(cor.kendall.ps_otu_adonis_ntf_nc.r.df_otu.ord.InterestOTUs)
kendall.OTUs.ps_otu.b.r = prune_taxa(kendall.OTUs, ps_otu.b.r)
tax_table(kendall.OTUs.ps_otu.b.r)
#The 5 OTUs are: n=3 Gallionella, n=1 Nitrospira, n=1 Omnitrophica (NA at Genus level)
# cors <- ddply(mtcars, c("vs", "am"), summarise, cor = round(cor(mpg, wt), 2))

#just for plotting purposes, change name of OTUs
# taxa_names(kendall.OTUs.ps_otu.b.r) <- paste0("OTU_", seq(ntaxa(kendall.OTUs.ps_otu.b.r)))
# kendall.OTUs.ps_otu.b.r.glom <- tax_glom(kendall.OTUs.ps_otu.b.r, taxrank = 'Species', NArm = FALSE)
kendall.OTUs.ps_otu.b.r.glom.ps_otudf_otu <- data.table(psmelt(kendall.OTUs.ps_otu.b.r))
kendall.OTUs.ps_otu.b.r.glom.ps_otudf_otu$Species <- as.character(kendall.OTUs.ps_otu.b.r.glom.ps_otudf_otu$Species)

kendall.OTUs.ps_otu.b.r.glom.ps_otudf_otu.plot<-ggplot(kendall.OTUs.ps_otu.b.r.glom.ps_otudf_otu[Abundance > 0], 
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

kendall.OTUs.ps_otu.b.r.glom.ps_otudf_otu.plot
# annotate("text", x=14, y=0.02, label = "R^2=0.78") +
# annotate("text", x=14, y=0.01, label = "alpha=0.00") +
# annotate("text", x=14, y=0.03, label = "beta=0.67")

ggsave("sgg_otu-KendallsCorSVs-TempOnly.plot.png", 
       height=8, 
       width=20)

# Prevalence in the dataset ----
#Completely copied from http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html#filter-taxa
mdt_otu = fast_melt(ps_otu.b)
prevdt_otu = mdt_otu[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = TaxaID]
#Test
# ggplot(prevdt_otu, aes(Prevalence)) + 
#   geom_histogram() + 
#   ggtitle("Histogram of Taxa Prevalence")
prevdt_otu[(Prevalence <= 0), .N]
prevdt_otu[(Prevalence <= 1), .N]
prevdt_otu[(Prevalence <= 2), .N]

# taxa cumulative sum
prevcumsum_otu = prevdt_otu[, .N, by = Prevalence]
setkey(prevcumsum_otu, Prevalence)
prevcumsum_otu[, CumSum := cumsum(N)]
pprevcumsum_otu = ggplot(prevcumsum_otu, aes(Prevalence, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Prevalence") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")
pprevcumsum_otu
#
addPhylum = unique(copy(mdt_otu[, list(TaxaID, Phylum)]))
# Join by TaxaID
setkey(prevdt_otu, TaxaID)
setkey(addPhylum, TaxaID)
prevdt_otu <- addPhylum[prevdt_otu]
showPhyla = prevdt_otu[, sum(TotalCounts), by = Phylum][order(-V1)][1:9]$Phylum
setkey(prevdt_otu, Phylum)
prev_totcounts<-ggplot(prevdt_otu[showPhyla], 
                       mapping = aes(Prevalence, TotalCounts, color = Phylum)) + 
  geom_point(size = 4, alpha = 0.7) + 
  scale_y_log10()
prev_totcounts
ggsave("sgg_otu-PrevalenceVSTotalCounts.plot.png", 
       height=8, 
       width=20)
# xxxx Does SILVA_132 influence general tax? YEAH it does ----
#waiting for DADA2 implementation