
#Loading required libraries
library(dada2)
library(decontam)


getwd()

# What do do first? modify paths, modify fastq_names, modify primers, 
# modify fragment lenght, configurate SILVA path in assign taxonomy
# 

#Setting paths
#path is the output_folder
#data_path is fastq_files_folder

path <- "C:/Users/danie/Desktop/MiSeq_SOP/results/"
data_path <- "C:/Users/danie/Desktop/MiSeq_SOP/"

#Getting fastq files and extracting sample names .... rename accordly with the sequences
fnFs <- sort(list.files(data_path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(data_path, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Checking data quality, filtering, and trimming .... 4 is the number of files to show quality of the reads


plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "CCGYCAATTYMTTTRAGTTT"
FWD_LENGTH <- nchar(FWD)
REV_LENGTH <- nchar(REV)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# The standard filtering parameters are starting points, not set in stone. 
# If you want to speed up downstream computation, consider tightening maxEE

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn the Error Rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
error_rate <- plotErrors(errF, nominalQ=TRUE)
error_rate <- plotErrors(errR, nominalQ=TRUE)



# Dereplication

derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- samples

# Inferring ASVs

dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo")
# dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo")
# dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=TRUE) # problem running this way if on Binder


# Merging Forward and Reverse reads

merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, trimOverhang=TRUE, minOverlap=170)

# this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
class(merged_amplicons) # list
length(merged_amplicons) # 20 elements in this list, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list 

class(merged_amplicons$B1) # each element of the list is a dataframe that can be accessed and manipulated like any ordinary dataframe

names(merged_amplicons$B1) # the names() function on a dataframe gives you the column names
# "sequence"  "abundance" "forward"   "reverse"   "nmatch"    "nmismatch" "nindel"    "prefer"    "accept"


# Generating a count table

seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) # 20 2521


# Chimera identification

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) # Identified 17 bimeras out of 2521 input sequences.

# though we only lost 17 sequences, we don't know if they held a lot in terms of abundance, this is one quick way to look at that
sum(seqtab.nochim)/sum(seqtab) # 0.9931372 # good, we barely lost any in terms of abundance


# Overview of counts throughout

# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

head(summary_tab)


# Assigning taxonomy

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
# creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)



# Removing likely contaminants ..... decontam package









































#Checking if primers were removed correclty by filterAndTrim
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
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs[[1]]), #or fn = fnFs.filtN
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtRs[[1]]), #or fn = fnRs.filtN
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs[[1]]), #or fn = fnFs.filtN
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtRs[[1]])) #or fn = fnRs.filtN

#Cutadapt is only needed when primers are not removed using filterAndTrim (e.g. variable length amplicons/primers or overlapping amplicons)
cutadapt <- "<path_to_cutadapt.exe_file>"
system2(cutadapt, args = "--version")
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                             fnFs.filtN[i], fnRs.filtN[i])) 
}
#Cutadapt over

#Leraning error rates
errF <- learnErrors(filtFs, nbases=1e+10, randomize=T, multithread=T)
errR <- learnErrors(filtRs, nbases=1e+10, randomize=T, multithread=T)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Dereplicating samples
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Run dada2 core algorithm
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[3]]
head(getUniques(dadaFs[[3]]))
head(getSequences(dadaFs[[3]]))
dadaFs_pool <- dada(derepFs, err=errF, pool = T, multithread=TRUE)
dadaRs_pool <- dada(derepRs, err=errR, pool = T, multithread=TRUE)
dadaFs_pool[[3]]
head(getUniques(dadaFs_pool[[3]]))
head(getSequences(dadaFs_pool[[3]]))

#Merging peired-end reads
mergers <- mergePairs(dadaFs_pool, derepFs, dadaRs_pool, derepRs, verbose=TRUE)
head(mergers[[3]])

#Creating ASV Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(311,416)]
dim(seqtab2)
table(nchar(getSequences(seqtab2)))

#Removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="pooled", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Tracking removed reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs_pool, getN), sapply(dadaRs_pool, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF_pooled", "denoisedR_pooled", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assigning taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "<path_to_silva_nr_v132_train_set.fa.gz_training_set_file>", multithread=TRUE)
taxa <- addSpecies(taxa, "<path_to_silva_species_assignment_v132.fa.gz_taxonomy_file>")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

#Converting dada2 output for use in Qiime2
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

#Creating fasta file of ASV sequences
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

#Make sure you open the generated .txt files and add the ASV_ID header in your ASV counts and taxa table