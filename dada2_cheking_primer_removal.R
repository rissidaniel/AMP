#Checking if primers were removed correclty by filterAndTrim

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
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs[[1]]), #or fn = fnFs.filtN
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtRs[[1]]), #or fn = fnRs.filtN
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs[[1]]), #or fn = fnFs.filtN
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtRs[[1]])) #or fn = fnRs.filtN



# If dada2 was not able to remove primers, use cutadapt

## Otherwise, contiue to Leraning error rates
