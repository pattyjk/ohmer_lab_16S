if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.20")
#resart or reload session

#load dada2
library(dada2)

#change to directory where fastq are (optional)
setwd("Documents/GitHub/ohmer_lab_16S/")

#get lsit of F/R reads
fnFs <- sort(list.files("./ohmer_demux/", pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files("./ohmer_demux/", pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#look at quality of reads (change the number to see different samples)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

#assign files names
filtFs <- file.path("./ohmer_demux/", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("./ohmer_demux/", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filter samples
filtered_out <- filterAndTrim(fnFs, filtFs,
                              fnRs, filtRs, 
                              rm.phix=TRUE, minLen=100, trimRight=10, trimLeft = 15)

#learn error rates on subset of data (use default)
errF <- learnErrors(filtFs, multithread=TRUE)
#100715374 total bases in 672839 reads from 73 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=TRUE)
#100288827 total bases in 672839 reads from 73 samples will be used for learning the error rates.

#optinal plot of error, uncomment
#plotErrors(errF, nominalQ=TRUE)

#do the dadaing for each read
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#merge pairs
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#make ASV table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#320 27261

#remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#Identified 18 bimeras out of 2141 input sequences.

#assign taxonomy
#download from: https://zenodo.org/records/4587955
taxa <- assignTaxonomy(seqtab.nochim, "./silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa2<-as.data.frame(taxa)

#make seq tab into R table
ohmer_lab_table<-as.data.frame(seqtab.nochim)

#write to file
write.table(ohmer_lab_table, 'asv_table.txt', quote=F, sep='\t')
write.table(taxa2, 'taxonomy.txt', quote=F, sep='\t')

#do a quick analysis
#load packages
library(vegan)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

#read in mapping file
meta<-read.delim('ohmer_lab_map.txt', header = T)

#set seed for reproducible analyses
set.seed(515)

#get sequencing depth
length(which(rowSums(ohmer_lab_table)<1000))
#4

#remove samples with <1000 reads
ohmer_lab_table<-ohmer_lab_table[-which(rowSums(ohmer_lab_table)<1000),]
dim(ohmer_lab_table)
#[1]   316 26233

#rarefy data to 1000 per sample
ohmer_rare<-rrarefy(ohmer_lab_table, sample=1000)

#calculate PCoA based on BC similarity
ohmer_pcoa<-capscale(ohmer_rare  ~ 1, distance='jaccard')

#pull out x/y coordinates
ohmer.scores<-scores(ohmer_pcoa)

#grab only sample coordinates, write to data frame
ohmer.coords<-as.data.frame(ohmer.scores$sites)

#create sample names as a column
ohmer.coords$SampleID<-row.names(ohmer.coords)

#map back meta data
ohmer.coords<-merge(ohmer.coords, meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained for first two axis
100*round(ohmer_pcoa$CA$eig[1]/sum(ohmer_pcoa$CA$eig), 3)
#10.5
100*round(ohmer_pcoa$CA$eig[2]/sum(ohmer_pcoa$CA$eig), 3)
#6.3

#plot PCoA
ggplot(ohmer.coords, aes(MDS1, MDS2, shape=Species, color=Type))+
  geom_point(size=2.8)+
  #geom_text()+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 10.5%")+
  ylab("PC2- 6.3%")

#calculate richness and Shannon
ohmer.alph<-as.data.frame(specnumber(ohmer_rare))
ohmer.alph$SampleID<-row.names(ohmer.alph)
ohmer.alph<-merge(ohmer.alph, meta, by='SampleID')

ohmer.alpha2<-as.data.frame(vegan::diversity(ohmer_rare, index = 'shannon'))
names(ohmer.alpha2)<-"Shannon"
ohmer.alph<-cbind(ohmer.alph, ohmer.alpha2)

#plot richness and shannon
ggplot(ohmer.alph, aes(Species, Shannon, fill=Type))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("Shannon Diversity")

ggplot(ohmer.alph, aes(Species, `specnumber(ohmer_rare)`, fill=Type))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("Shannon Diversity")
