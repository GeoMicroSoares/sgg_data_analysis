#completely adapted from https://benjjneb.github.io/dada2/tutorial.html
library(dada2); packageVersion("dada2")
path <- "/media/andre/B2F8C9A0F8C962E9/SGG_16S_analysis/SGG_16S_data"
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# Place filtered files in filtered/ subdirectory
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
#filter and trim reads
#F reads were truncated at 260bp and R at 200 based on:
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
#
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(260, 200), maxN=0, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
# On Windows set multithread=FALSE
head(out)
#learn errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#merger with F and R reads didn't work. alternative with F reads only:
filts <- list.files(filt_path, pattern="fastq.gz", full.names=TRUE) # CHANGE if different file extensions
sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1) # Assumes filename = sample_XXX.fastq.gz
names(filts) <- sample.names
set.seed(100)
err <- learnErrors(filts, nreads = 1e6, multithread=TRUE, randomize=TRUE)
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}
# Construct sequence table and write to disk
seqtab <- makeSequenceTable(dds)

###

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#sanity check for seq lenght
table(nchar(getSequences(seqtab)))
#this could help to getter a better handle of their sizes
plot(table(nchar(getSequences(seqtab))))
#remove bimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)