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

seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/seqtab.rds")
# Remove chimeras
seqtab_bimR <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
saveRDS(seqtab_bimR, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/seqtab_final.rds")
# seqtab_bimR<-readRDS("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/seqtab_final.rds")
# Assign taxonomy, add species
tax <- assignTaxonomy(seqtab_bimR, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SILVA_128/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(tax, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SILVA_128/silva_species_assignment_v128.fa.gz")
saveRDS(taxa, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/tax_final.rds")
#test
taxa.print <- taxa
rownames(taxa.print) <- NULL
#see if tax assignments worked
head(taxa.print)

#Check SV table size with phyloseq
#taxa_are_rows bc of makeSequenceTable()
sgg_metadata<-read.csv("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_metadata/SGG_metadata.csv", header=TRUE)
#reorder sgg_metadata based on seqtab_bimR
order<-row.names(seqtab_bimR)
sgg_metadata_ord<-sgg_metadata[match(order, sgg_metadata$Sample_Code),]
rownames(sgg_metadata_ord) <- sgg_metadata_ord$Sample_Code

ps <- phyloseq(otu_table(seqtab_bimR, taxa_are_rows=FALSE),
               sample_data(sgg_metadata_ord), 
               tax_table(taxa))

ps

# 0 - Basic stats 
#DADA2 processing yields
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab_bimR))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)

#In particular: how many reads left in seqtab_bimR?
sum(rowSums(seqtab_bimR))
#How many OTUs in ps?
ntaxa(ps)
#70,937 SVs to 9,656,492 reads

#Seq Depth
psdt = data.table(as(sample_data(ps), "data.frame"),
                 TotalReads = sample_sums(ps), keep.rownames = TRUE)
setnames(psdt, "rn", "SampleID")
pSeqDepth = ggplot(psdt, aes(TotalReads)) + 
  geom_histogram(binwidth = 4000) + 
  ggtitle("Sequencing Depth")
pSeqDepth

ggsave("pSeqDepth.png", width=14, height=6)

#Number of SVs per phylum
sv.phyl<-as.data.frame(table(tax_table(ps)[, "Phylum"], exclude = NULL))
sv.phyl.ord <- sv.phyl[order(-sv.phyl$Freq),] 
# uh la la Parcubacteria!

#Number of sequences per SV, Sequences per sample
readsumsdf = data.frame(nreads = sort(taxa_sums(ps), TRUE), sorted = 1:ntaxa(ps), 
                        type = "SVs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps), 
                                                        TRUE), sorted = 1:nsamples(ps), type = "Samples"))
p_read_sums = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p_read_sums + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
ggsave("SVsSeq_SeqsSpl.png", width=14, height=6)

#General phyla table:
#transformed dataset
ps.phylum <- ps
ps.phyl.glom <- tax_glom(ps.phylum,taxrank = "Phylum") 
sample_data(ps.phyl.glom)$Month = factor(sample_data(ps.phyl.glom)$Month, 
                                         levels = c("April","August","December","Control"))
ps.phyl.glom.a = subset_samples(ps.phyl.glom, Site.name != "Taffs Well PUMPED")
ps.phyl.glom.b = subset_samples(ps.phyl.glom.a, Sample_type != "Control")

ps.phyl.glom.r <-  transform_sample_counts(ps.phyl.glom.b, function(x) {x/sum(x)} ) 
ps.phyl.glom.r.m <- psmelt(ps.phyl.glom.r)       

ps.phyl.glom.r.m.a.ord.plot<-ggplot(ps.phyl.glom.r.m.a.ord, aes(x = Site.name, y = Abundance/3, fill = Phylum)) + 
  facet_grid(Month~.) +
  geom_bar(stat = "identity") +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  ylab("Relative Abundance\n")
ggsave("ps.phyl.glom.r.m.a.ord.allphylatrans.png", height=8, width=16)

##all phyla
##top 10 phyla, rest = "Other"
##transformed dataset
ps.phyl.glom.b.r <-  transform_sample_counts(ps.phyl.glom.b, function(x) {x/sum(x)} ) 
ps.phyl.glom.ba <- psmelt(ps.phyl.glom.b.r)       
ps.phyl.glom.ba$Phylum <- as.character(ps.phyl.glom.ba$Phylum)
# ps.phyl.glom.ba$Abundance <- as.(ps.phyl.glom.ba$Abundance)
# group dataframe by Phylum, calculate median rel. abundance
medians <- ddply(ps.phyl.glom.ba, ~Phylum, function(x) c(median=median(x$Abundance)))
#find phyla less abundant than the top10
other_top10<-order(medians$median, decreasing = TRUE)
other_top10<-medians[order(-medians$median), , drop = FALSE]
other_top10.c<-tail(other_top10,-10)
# change their name to "Other"
ps.phyl.glom.ba[ps.phyl.glom.ba$Phylum %in% other_top10.c$Phylum,]$Phylum <- 'Other'

cols_sgg <- colorRampPalette(brewer.pal(8, "Set1"))
myPal_sgg <- cols(length(unique(ps.phyl.glom.ba$Phylum)))

ps.phyl.glom.plot<-ggplot(ps.phyl.glom.ba, aes(x = Site.name, y = Abundance/3, fill = Phylum)) + 
  facet_grid(Month~.) +
  geom_bar(stat = "identity") +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  # coord_cartesian(ylim = c(0, 1)) +
  ylab("Relative Abundance\n")
  # scale_color_manual(values=myPal_sgg)
ps.phyl.glom.ba.plot
ggsave("ps.phyl.glom.r.m.a.ord.allphylatrans.top10.png", height=8, width=16)

##all phyla
##top 10 phyla, rest = "Other"
##NOT transformed dataset
ps.phyl.glom.bp <- psmelt(ps.phyl.glom.b)       
ps.phyl.glom.bp$Phylum <- as.character(ps.phyl.glom.bp$Phylum)
# ps.phyl.glom.ba$Abundance <- as.(ps.phyl.glom.ba$Abundance)
# group dataframe by Phylum, calculate median rel. abundance
medians.bp <- ddply(ps.phyl.glom.bp, ~Phylum, function(x) c(median=median(x$Abundance)))
#find phyla less abundant than the top10
other_top10.bp<-order(medians.bp$median, decreasing = TRUE)
other_top10.bp<-medians.bp[order(-medians.bp$median), , drop = FALSE]
other_top10.bp.c<-tail(other_top10.bp,-10)
# change their name to "Other"
ps.phyl.glom.bp[ps.phyl.glom.bp$Phylum %in% other_top10.bp.c$Phylum,]$Phylum <- 'Other'

ps.phyl.glom.bp.p<-ps.phyl.glom.bp

ps.phyl.glom.bp.p.plot<-ggplot(ps.phyl.glom.bp.p, aes(x = Site.name, y = Abundance/3, fill = Phylum)) + 
  facet_grid(Month~.) +
  geom_bar(stat = "identity") +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  # coord_cartesian(ylim = c(0, 1)) +
  ylab("Relative Abundance\n")
# scale_color_manual(values=myPal_sgg)
ps.phyl.glom.bp.p.plot
ggsave("ps.phyl.glom.r.m.a.ord.allphylanottrans.top10.png", height=8, width=16)

##all phyla
##top 10 phyla, rest = "Other"
##transformed dataset
##proteobacterialclasses?
ps.prot<-ps
sample_data(ps.prot)$Month = factor(sample_data(ps.prot)$Month, 
                                         levels = c("April","August","December","Control"))
ps.prot.a = subset_samples(ps.prot, Site.name != "Taffs Well PUMPED")
ps.prot.b = subset_samples(ps.prot.a, Sample_type != "Control")

ps.phyl.glom.b.r <-  transform_sample_counts(ps.prot.b, function(x) {x/sum(x)} )
ps.phyl.glom.r.prot = subset_taxa(ps.phyl.glom.b.r, Phylum == "Proteobacteria")
ps.phyl.glom.r.prot.cl <- tax_glom(ps.phyl.glom.r.prot, taxrank="Class")
# ps.phyl.glom.r.prot.cl.t <-  transform_sample_counts(ps.phyl.glom.r.prot.cl, function(x) {x/sum(x)} )
ps.phyl.glom.r.prot.cl.p <- psmelt(ps.phyl.glom.r.prot.cl)
ps.phyl.glom.r.prot.cl.p$Class <- as.character(ps.phyl.glom.r.prot.cl.p$Class)
# group dataframe by Phylum, calculate median rel. abundance
medians.cl.p <- ddply(ps.phyl.glom.r.prot.cl.p, ~Class, function(x) c(medians=median(x$Abundance)))
#find classes less abundant than the top10
other_top10.cl<-medians.cl.p[order(-medians.cl.p$median), , drop = FALSE]
other_top10.cl.c<-tail(other_top10.cl,-5)
# change their name to "Other"
ps.phyl.glom.r.prot.cl.p[ps.phyl.glom.r.prot.cl.p$Class %in% other_top10.cl.c$Class,]$Class <- 'Other Proteobacteria'
#remove proteobacteria from this

other_top5.c<-tail(other_top10.bp,-15)
# change their name to "Other"
ps.phyl.glom.ba.top5<-ps.phyl.glom.ba
ps.phyl.glom.ba.top5[ps.phyl.glom.ba.top5$Phylum %in% other_top5.c$Phylum,]$Phylum <- 'Other'

ps.phyl.glom.ba.top5.noP<-ps.phyl.glom.ba.top5[!grepl("Proteobacteria", ps.phyl.glom.ba.top5$Phylum),]
#paste these two together but change Class or Phylum to Taxa
names(ps.phyl.glom.r.prot.cl.p)[names(ps.phyl.glom.r.prot.cl.p) == 'Class'] <- 'Taxa'
#remove Phylum
ps.phyl.glom.r.prot.cl.p.nP<-ps.phyl.glom.r.prot.cl.p[ , -which(names(ps.phyl.glom.r.prot.cl.p) %in% c("Phylum"))]
names(ps.phyl.glom.ba.top5.noP)[names(ps.phyl.glom.ba.top5.noP) == 'Phylum'] <- 'Taxa'
ps.phyl.glom.ba.top5.noP.trim<-ps.phyl.glom.ba.top5.noP[ , -which(names(ps.phyl.glom.ba.top5.noP) %in% c("Class","Order",
                                                                                                         "Family","Genus",
                                                                                                         "Species"))]
ps.phy.prot.cl <- rbind(ps.phyl.glom.r.prot.cl.p.nP, ps.phyl.glom.ba.top5.noP.trim)
 
ps.phy.prot.cl.plot<-ggplot(ps.phy.prot.cl, aes(x = Site.name, y = Abundance/3, fill = Taxa)) + 
  facet_grid(Month~.) +
  geom_bar(stat = "identity") +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  # coord_cartesian(ylim = c(0, 1)) +
  ylab("Relative Abundance\n")
ps.phy.prot.cl.plot

ggsave("ps.phyl.glom.r.m.a.ord.allphylanottrans.png", height=8, width=16)


#all phyla
#not transformed dataset
ps.phyl.glom.b.r.m <- psmelt(ps.phyl.glom.b)       
ps.phyl.glom.b.r.m.plot<-ggplot(ps.phyl.glom.b.r.m, aes(x = Site.name, y = Abundance/3, fill = Phylum)) + 
  facet_grid(Month~.) +
  geom_bar(stat = "identity") +
  # Remove x axis title
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  # coord_cartesian(ylim = c(0, 1)) +
  ylab("Relative Abundance\n")
ggsave("ps.phyl.glom.r.m.a.ord.allphylanottrans.png", height=8, width=16)

# 1 - Tax profiles across month and sites
#Proteobacteria barplots
#subset only Proteobacteria
ps.proteob = subset_taxa(ps, Phylum=="Proteobacteria")
#here new order is given to some variables, for plotting reasons
sample_data(ps.proteob)$Month = factor(sample_data(ps.proteob)$Month, 
                                       levels = c("April","August","December","Control"))
ps.proteob.noTFP = subset_samples(ps.proteob, Site.name != "Taffs Well PUMPED")
ps.proteob.noTFP.noCtrl = subset_samples(ps.proteob.noTFP, Sample_type != "Control")

ps.proteob.noTFP.noCtrl.bar<-plot_bar(ps.proteob.noTFP.noCtrl, 
                                      fill="Class", facet_grid=~Month, x="Site.name") + 
  geom_bar(aes(color=Class, fill=Class), 
           stat="identity", position="stack")
ps.proteob.noTFP.noCtrl.bar
ggsave("ps.proteob.noTFP.noCtrl.bar.png",width = 18,height=8)

#Proteobacterial orders barplots - beta and epsilonproteobacteria
#subset only Proteobacteria
ps.proteob.cl = subset_taxa(ps, Class=="Betaproteobacteria" | Class=="Epsilonproteobacteria")
#here new order is given to some variables, for plotting reasons
sample_data(ps.proteob.cl)$Month = factor(sample_data(ps.proteob.cl)$Month, 
                                          levels = c("April","August","December","Control"))
ps.proteob.cl.noTFP = subset_samples(ps.proteob.cl, Site.name != "Taffs Well PUMPED")
ps.proteob.cl.noTFP.noCtrl = subset_samples(ps.proteob.cl.noTFP, Sample_type != "Control")

ps.proteob.cl.noTFP.noCtrl.bar<-plot_bar(ps.proteob.cl.noTFP.noCtrl, 
                                         fill="Order", facet_grid=~Month, x="Site.name") + 
  geom_bar(aes(color=Order, fill=Order), 
           stat="identity", position="stack")
ps.proteob.cl.noTFP.noCtrl.bar
ggsave("ps.proteob.cl.noTFP.noCtrl.bar.png",width = 18,height=8)

#Proteobacterial orders barplots - Gallionellaceae family
ps.proteob.g = subset_taxa(ps, Family=="Gallionellaceae")

sample_data(ps.proteob.g)$Month = factor(sample_data(ps.proteob.cl)$Month, 
                                          levels = c("April","August","December","Control"))
ps.proteob.g.noTFP = subset_samples(ps.proteob.g, Site.name != "Taffs Well PUMPED")
ps.proteob.g.noTFP.noCtrl = subset_samples(ps.proteob.g.noTFP, Sample_type != "Control")

ps.proteob.g.noTFP.noCtrl.bar<-plot_bar(ps.proteob.g.noTFP.noCtrl, 
                                         fill="Genus", facet_grid=~Month, x="Site.name") + 
  geom_bar(aes(color=Genus, fill=Genus), 
           stat="identity", position="stack")
ps.proteob.g.noTFP.noCtrl.bar
ggsave("ps.proteob.g.noTFP.noCtrl.bar.png",width = 18,height=8)

#Proteobacterial orders barplots - Campylobacterales order
ps.proteob.c = subset_taxa(ps, Order=="Campylobacterales")

sample_data(ps.proteob.c)$Month = factor(sample_data(ps.proteob.cl)$Month, 
                                          levels = c("April","August","December","Control"))
ps.proteob.c.noTFP = subset_samples(ps.proteob.c, Site.name != "Taffs Well PUMPED")
ps.proteob.c.noTFP.noCtrl = subset_samples(ps.proteob.c.noTFP, Sample_type != "Control")

ps.proteob.c.noTFP.noCtrl.bar<-plot_bar(ps.proteob.c.noTFP.noCtrl, 
                                         fill="Family", facet_grid=~Month, x="Site.name") + 
  geom_bar(aes(color=Family, fill=Family), 
           stat="identity", position="stack")
ps.proteob.c.noTFP.noCtrl.bar
ggsave("ps.proteob.c.noTFP.noCtrl.bar.png",width = 18,height=8)

# 2 - nMDS on dataset
# #rename SVs to make it easier for nMDS to take it?
# taxa_names(ps.ord) <- paste0("SV", seq(ntaxa(ps.ord)))
#empty control samples were impeding calculations by giving NAs like crazy. out now!
ps.nc = subset_samples(ps, Sample_type != "Control")
ps.nc.r = transform_sample_counts(ps.nc, function(x) x/sum(x))

dist_methods <- unlist(distanceMethodList)
print(dist_methods)
dist_methods <- dist_methods[-(1:3)]
dist_methods["designdist"]
dist_methods = dist_methods[-which(dist_methods=="ANY")]

plist <- vector("list", length(dist_methods))
names(plist) = dist_methods

for( i in dist_methods ){
  # Calculate distance matrix
  ps.nc.r.dist <- phyloseq::distance(ps.nc.r, method=i)
  # Calculate ordination
  ps.nc.r.nmds  <- ordinate(ps.nc.r, "MDS", distance=ps.nc.r.dist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(ps.nc.r, ps.nc.r.nmds, color="SeqTech", shape="Enterotype")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}

df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Site.name, shape=Month))
p = p + geom_point(size=3, alpha=0.5)
p = p + facet_wrap(~distance, scales="free")
p = p + ggtitle("MDS on various distance metrics for SGG dataset")
p

plot_scree(ps.nc.r.nmds, "Scree Plot: Bray-Curtis MDS")

#3 - Does temperature relate to any of the chars of the microbiome?
#a) plot temperature across sites and month
ps_temp<-ps
sample_data(ps_temp)$Month = factor(sample_data(ps_temp)$Month, 
                                          levels = c("April","August","December","Control"))
sample_data(ps_temp)$Site.name = factor(sample_data(ps_temp)$Site.name, 
                                    levels = c("Blaenavon","Cefn Hengoed","Mountain Gate","Taff Bargoed",
                                               "Glyncastle","Dinas","Ynysarwed","Celynen North","Morlais","Lindsay",
                                               "Six Bells","Crumlin Navigation","Taffs Well"))
ps_temp.tw = subset_samples(ps_temp, Site.name != "Taffs Well PUMPED")
ps_temp.tw.c = subset_samples(ps_temp.tw, Sample_type != "Control")
temp_table<-sample_data(ps_temp.tw.c)[,c("Sample_Code","Month","Site.name","Temp.C")]
ggplot(temp_table, aes(x=Site.name, y=Temp.C, fill=Month)) +
  geom_boxplot() +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1)) +
  geom_hline(yintercept = 13.4,
             color="dark red",
             linetype="dashed") +
  annotate("text", "Mountain Gate", 13.5, vjust = -1, 
           label = "Average mine water temperature (Farr et al, 2016)") +
  labs(y=expression("Temperature in"*~degree*C)) +
  geom_ribbon(aes(ymin=9, ymax=11), fill = "grey70") +
  ggtitle("Temperatures across the SGG sampled sites")
ggsave("Temperatures across the SGG sampled sites.png")

#b) relating OTUs, clades to temp?
#4 - adonis to relate SVs to env variables
sgg_metadata_adonis<-read.csv("/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_metadata/SGG_metadata_adonis.csv", header=TRUE)
#reorder sgg_metadata based on seqtab_bimR
order<-row.names(seqtab_bimR)
sgg_metadata_adonis_ord<-sgg_metadata_adonis[match(order, sgg_metadata_adonis$Sample_Code),]
rownames(sgg_metadata_adonis_ord) <- sgg_metadata_adonis_ord$Sample_Code

ps_adonis <- phyloseq(otu_table(seqtab_bimR, taxa_are_rows=FALSE),
               sample_data(sgg_metadata_adonis_ord), 
               tax_table(taxa))
ps_adonis

sample_data(ps_adonis)$Month = factor(sample_data(ps_adonis)$Month, 
                                         levels = c("April","August","December","Control"))
ps_adonis_ntf = subset_samples(ps_adonis, Site.name != "Taffs Well PUMPED")
ps_adonis_ntf_nc = subset_samples(ps_adonis_ntf, Sample_type != "Control")

df = as(sample_data(ps_adonis_ntf_nc), "data.frame")
d = phyloseq::distance(ps_adonis_ntf_nc, "bray")
ps.adonis.ntf.nc.adonis = adonis(d ~ Month + Site.name + Temp.C + pH + EC.uS.cm, df)

ps.adonis.ntf.nc.adonis
plot(ps.adonis.ntf.nc.adonis$aov.tab)
densityplot(permustats(ps.adonis.ntf.nc.adonis))

#now from https://rpubs.com/dillmcfarlan/R_microbiotaSOP
#see what SVs correlate to temperature with cor()
ps_adonis_ntf_nc.r = transform_sample_counts(ps_adonis_ntf_nc, function(x) x/sum(x))
SV1 = as(otu_table(ps_adonis_ntf_nc.r), "matrix")
# transpose if necessary
if(taxa_are_rows(ps_adonis_ntf_nc.r)){SV1 <- t(SV1)}
# Coerce to data.frame
SV1df = as.data.frame(SV1)

cor.spearman.ps_adonis_ntf_nc.r = cor(SV1df, 
                                      sample_data(ps_adonis_ntf_nc.r)$Temp.C, 
                                      method = "spearman")

cor.kendall.ps_adonis_ntf_nc.r = cor(SV1df, 
                                      sample_data(ps_adonis_ntf_nc.r)$Temp.C, 
                                      method = "kendall")

layout(matrix(c(1,1), 1, 1, byrow = TRUE))
plot(sort(cor.spearman.ps_adonis_ntf_nc.r, decreasing = TRUE))
abline(h = c(-0.5,0, 0.5))
plot(sort(cor.kendall.ps_adonis_ntf_nc.r, decreasing = TRUE))
abline(h = c(-0.5,0, 0.5))

#look into this!
beta <- betadisper(d, df$Temp.C)
permutest(beta)

#5 - network
ps.network <- ps
sample_data(ps.network)$Month = factor(sample_data(ps.network)$Month, 
                                         levels = c("April","August","December","Control"))
ps.network.a = subset_samples(ps.network, Site.name != "Taffs Well PUMPED")
ps.network.b = subset_samples(ps.network.a, Sample_type != "Control")
ps.network.b.r = transform_sample_counts(ps.network.b, function(x) x/sum(x))

# ig <- make_network(ps.network.b.r, dist.fun="bray", max.dist = 0.4)
sgg_net<-plot_net(ps.network.b, distance = "bray", maxdist = 0.7, laymeth = "fruchterman.reingold",
         color="Site.name", shape="Month",
         point_alpha=0.9, point_size=4) 
sgg_net +
  ggtitle("Bray-Curtis distance Fruchterman-Reingold network")
