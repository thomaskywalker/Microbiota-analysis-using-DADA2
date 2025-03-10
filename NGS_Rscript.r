###if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.14")
#install.packages("Rcpp")
##

#library packages
library("Rcpp")
library("dada2")
library(dada2); packageVersion("dada2")
library("dada2")
library("stats")
library("ggplot2")

#set file path 
path <- "C:/Users/user/Desktop/NGS"
list.files(path)

#check file content
library(ShortRead)
readFastq("C:/Users/user/Desktop/NGS/010SRR1260406.fastq")
readFastq("C:/Users/user/Desktop/NGS/009SRR1260406.fastq")

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern=".fastq"))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
plotQualityProfile(fnFs[1:3])
??plotQualityProfile

#Perform filtering and trimming
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, ".fastq"))
out <- filterAndTrim(fnFs, filtFs,trimRight=15, trimLeft = 15,minLen = 20,  maxN=0, truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=F)
?filterAndTrim
head(out)
out

plotQualityProfile(filtFs[7:9])

# Start the clock
ptm <- proc.time()
# Learn error rates
errF <- learnErrors(filtFs, nbases = 1e+10, multithread=F)
#stop the clock
proc.time() - ptm
plotErrors(errF, nominalQ=TRUE)

#Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

#sample inference
# Start the clock
ptm <- proc.time()
dadaFs <- dada(derepFs, err=errF, multithread=FALSE,HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
# Stop the clock
proc.time() - ptm
dadaFs[[1]]

#make sequence table directly out of dadaFs
#do the rest of analysis with seqtab
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
?dim

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)), main="Distribution of sequence lengths")

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=F, verbose=TRUE)
?removeBimeraDenovo
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
?cbind
colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
track

#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/user/Desktop/NGS/silva_nr_v132_train_set.fa.gz", multithread=F)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#export ASV table
setwd(path)
seqtable.taxa.plus <- cbind("seq"= rownames(taxa),t(seqtab.nochim),taxa)
write.table(seqtab.nochim,"dada2_counts.txt",sep="\t" , quote = F , row.names = T)
write.table(seqtable.taxa.plus,"dada2_counts.taxon.species.txt" , sep = "\t" , quote=F , row.names = F)
write.table(track,"dada2_track.txt" , sep = "\t" , quote=F , row.names = F)
write.csv(seqtab.nochim,"ASV.csv",sep=",",row.names=T, na = "NA")
write.csv(seqtable.taxa.plus,"dada2_counts.taxon.species.csv" , sep = "," , row.names = T)

#Evaluating DADA2’s accuracy on the mock community:
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

#BiocManager::install("phyloseq")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
# Make a data.frame holding the sample data
otutable=read.table("C:/Users/user/Desktop/NGS/dada2_counts.taxon.species.txt")

samples.out <- rownames(seqtab.nochim)
samples.out
subject <- sapply(strsplit(samples.out, "S"), `[`, 1)
subject
gender <- substr(subject,1,5)
gender
subject <- substr(subject,2,999)
substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Sample=subject)
rownames(samdf) <- samples.out

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
#ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
dna <- Biostrings::DNAStringSet(taxa_names(ps))
dna
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps <- subset_taxa(ps, !is.na(Family) & !Family %in% c("", "NA"))
ps

#Plot the species richness (alpha diversity)
dev.off()
plot_richness(ps, x="Sample", measures=c("Observed","Chao1","Shannon", "Simpson"), color="Sample") + theme_bw()
sample_data 

#Beta diversity PcoA
ord = ordinate(ps, method="PCoA", distance = "bray")
?plot_ordination
plot_ordination(ps.prop, ord,color = "sample_Sample",title = "Beta diversity") + 
  geom_point(size=3) 

#Create ordination plots
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps, ord.nmds.bray,  title="Bray NMDS")

#Bar plot
top50 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)
plot_bar(ps.top50, fill="Family")

#heatmap
sampleOrder = unique(sample_names(ps))
taxaOrder = unique(taxa_names(ps))
plot_heatmap(ps.top50, taxa.label="Family",sample.order=sampleOrder)
?plot_heatmap
