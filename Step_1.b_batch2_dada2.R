# Processing of raw reads, 2nd batch of 16S rRNA amplicons.
library(dada2); version()
setwd("batch2_data")

path<-"16S_amplicons"
list.files(path)
library(ShortRead)
library(Biostrings)
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[2:4], collapse = "_"))
plotQualityProfile(fnFs[3]) # looks okay until 280 nt
plotQualityProfile(fnRs[2]) #looks okay until 180 nt
table(nchar(getSequences(fnFs[2])))
#Searching the reads for primers (V5-V7)
FWD<-"AACMGGATTAGATACCCKG"
REV<-"ACGTCATCCCCACCTTCC"
nchar(FWD)
nchar(REV)
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
#create paths for filtered reads
filtFsNS <- file.path(path, "filtN", paste0(sample.names, "_F_filtN.fastq.gz"))
filtRsNS <- file.path(path, "filtN", paste0(sample.names, "_R_filtN.fastq.gz"))

names(filtFsNS) <- sample.names
names(filtRsNS) <- sample.names
#Trimming the reads to remove N's
out<-filterAndTrim(fnFs, filtFsNS, fnRs, filtRsNS, maxN = 0, truncLen = c(280,200),trimLeft =c(19,18), 
              rm.phix = TRUE, truncQ = 2, maxEE = c(2,2), compress = TRUE, multithread = TRUE)
#Large proportion of reads are kept! Looks good

#Here I am only checking the primers if they are found in reverse. Trucanting the reads, will remove them in opposite direction.
#Check length distribution
table(nchar(getSequences(filtFsNS[1])))
#Function for search of primers
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

out
errF <- learnErrors(filtFsNS, multithread=TRUE)
plotErrors(errF, nominalQ = TRUE) #error rates look okay
errR <- learnErrors(filtRsNS, multithread=TRUE) #
plotErrors(errR, nominalQ =TRUE)

dadaFs <-dada(filtFsNS, err = errF, multithread = TRUE) 
dadaRs<-dada(filtRsNS, err = errR, multithread = TRUE)

mergers <- mergePairs(dadaFs, filtFsNS, dadaRs, filtRsNS, verbose=TRUE)
head(mergers)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab))) #They are between 357 and 404 in length. I will keep all, and problematic ones will most 
#likely be removed as chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab,  multithread=TRUE, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab2)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim)/out[,1])
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim","ratio")
track

setwd("../cintia");taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE);setwd("../Frederik")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print[1:20,]

#at the dada2 tutorial there is a suggestion on how to look at the mock community specifically.

save.image("batch2_syncom.rData")
saveRDS(taxa,"taxa_batch2.rds")
saveRDS(seqtab.nochim,"asv_table.rds")

#Convert to a phyloseq object (easy to work with)
samples.out <- rownames(seqtab.nochim)

library(phyloseq)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
ps
dna <- Biostrings::DNAStringSet(taxa_names(ps)) #Here the DNA sequences are stored in a separate object, which is then put into the phyloseq object. It might be usefull later on in an anlysis.
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) #Change taxa names to ASV1, ASV2, etc.
ps
sample_names(ps)<-sample.names
