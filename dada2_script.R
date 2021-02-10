library(dada2); packageVersion("dada2")
library(phyloseq)

# MODIFY: path of fastq files and reference database; fastq name format if necessary; quality control parameters,
# primers in the primer check; metadata file (qiime1 format in the example); 

# to check metadata file https://keemei.qiime2.org/

# setwd("c:/Users/danie/Desktop/dada2_big_data/16S_chemerin_tutorial/fastq/")
path <- "c:/Users/danie/Desktop/dada2_big_data/16S_chemerin_tutorial/fastq/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


###Filter and trim

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,210),
                     maxN=0, maxEE=c(3,7), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
#head(out)

### checking if primers were correctly removed 
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "CCGYCAATTYMTTTRAGTTT"
FWD_LENGTH <- nchar(FWD)
REV_LENGTH <- nchar(REV)


allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
library("ShortRead")
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
primer_check <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs[[1]]), #or fn = fnFs.filtN
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtRs[[1]]), #or fn = fnRs.filtN
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs[[1]]), #or fn = fnFs.filtN
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtRs[[1]])) #or fn = fnRs.filtN

write.csv(primer_check,"primer_check.csv")



### Learn the Error Rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)

### Sample Inference

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

### Merge paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
# head(mergers[[1]])

### Construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


### Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

### Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "tracking_reads.csv")

### Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "c:/Users/danie/Downloads/silva_nr_v138_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "c:/Users/danie/Downloads/silva_species_assignment_v138.fa.gz")

#taxa.print <- taxa # Removing sequence rownames for display only
#rownames(taxa.print) <- NULL
#head(taxa.print)


######################################################################## create phyloseq object

library(phyloseq)
library(Biostrings)
library(ggplot2)

### saving dada2 objects for phyloseq scritps

saveRDS(seqtab.nochim, "seqtab.nochim.rds")
saveRDS(taxa, "taxa.rds")
save.image("my_session.rdata")

### to open object in new R script, in the case seqtab.nochim is the otu_table
seqtab.nochim <- readRDS("seqtab.nochim.rds")
taxa <- readRDS("taxa.rds")


samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)

samdf <- read.csv("map.txt", sep = "\t", row.names = 1)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))



dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


############################################################### Converting dada2 output for use in Qiime2


asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Creating fasta file of ASV sequences
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

#Creating a count table
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.txt", sep="\t", quote=F)

#Creating a taxonomy table
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.txt", sep="\t", quote=F)


