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

# Filename parsing
path <- "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed"
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
fns <- list.files(path, pattern="fastq.gz") # CHANGE if different file extensions

#check qual profiles before deciding trimming, etc
fnFs <- sort(list.files(path, pattern="_R1_001_trim.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001_trim.fastq", full.names = TRUE))

plotQualityProfile(fnFs[1:2]) + 
  scale_x_continuous(breaks=c(0,50,100,150,200,250,300)) +
  scale_y_continuous(breaks=c(0,10,20,22,24,26,28,30,35,40)) +
  theme(panel.grid.major = element_line(colour = "black"))
plotQualityProfile(fnRs[1:2])+ 
  scale_x_continuous(breaks=c(0,50,100,150,200,250,300)) +
  scale_y_continuous(breaks=c(0,10,20,22,24,26,28,30,35,40)) +
  theme(panel.grid.major = element_line(colour = "black"))

filt_path <- file.path(path, "filtered")
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt_trim.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt_trim.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(260, 200), maxN=0, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
#sample, sequence variant inference
filts <- list.files(filtpath, pattern="fastq.gz", full.names=TRUE) # CHANGE if different file extensions
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = sample_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = sample_XXX.fastq.gz
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# Learn error rates
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nread=1e6, multithread=TRUE, MAX_CONSIST = 30)
# Learn reverse error rates
errR <- learnErrors(filtRs, nread=1e6, multithread=TRUE, MAX_CONSIST = 20)
# Infer sequence variants
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Construct sequence table and write to disk
#"A row for each sample, and a column for each unique sequence across all the samples"
seqtab <- makeSequenceTable(mergers)
#lengths
table(nchar(getSequences(seqtab)))
saveRDS(seqtab, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/seqtab.rds")
# Remove chimeras
seqtab_bimR <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
saveRDS(seqtab_bimR, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/seqtab_final.rds")
# Assign taxonomy, add species
tax <- assignTaxonomy(seqtab_bimR, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SILVA_128/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(tax, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SILVA_128/silva_species_assignment_v128.fa.gz")
# Write to disk
saveRDS(taxa, "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data_trimmed_DADA2_table/tax_final.rds")
# #test
# taxa.print <- taxa
# rownames(taxa.print) <- NULL
# #see if tax assignments worked
# head(taxa.print)
#
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
