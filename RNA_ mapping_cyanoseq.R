# https://benjjneb.github.io/dada2/tutorial.html

library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions

### Load Libraries ####
library(data.table)
library(dada2)
library(beepr)

### Set Path ####
path <- "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/BSC Research/Ch. 4 Biocrust Burn/RNA/00_fastq"
list.files(path)

### Forward and reverse gz files have the format samplename_1 and samplename_2 ####
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))


### Extract sample names, assuming filenames have format: SAMPLENAME_XXX.gz####
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list.files(path)

###plot the read quality, these look good! Q score is above 30####
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

### Filter and Trim ####
# Place filtered files in filtered/ subdirectory
filtFs <- file.path("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/BSC Research/Ch. 4 Biocrust Burn/RNA/00_fastq/filtered", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/BSC Research/Ch. 4 Biocrust Burn/RNA/00_fastq/filtered", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
list.files(path)




# We’ll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.

# use the trim left function to remove the primers 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft =c(16,12)) #beep(sound = 3) # On Windows set multithread=FALSE
head(out)

### dereplication -- this reduces the downstrean computationsl time.  ####
# when I first ran this code, calculating the error rates took 12 hours so hopefully adding the dereplication step will speed up the process
# increase memory first 
derepF1 <- derepFastq(filtFs, verbose=TRUE); beep(sound = 3)  
derepR1 <- derepFastq(filtRs, verbose=TRUE); beep(sound = 3) 

### Learn the error rates ####

errF <- learnErrors(filtFs, multithread = TRUE)#beep(sound = 3) # multithread does not work on windows :(
errR <- learnErrors(filtRs, multithread = TRUE)# beep(sound = 3)
plotErrors(errF, nominalQ=TRUE)

### Sample inference ####
dadaFs <- dada(filtFs, err=errF, multithread = TRUE)# beep(sound = 3)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)# beep(sound = 3)
dadaFs[[1]]

### merged paired reads ####

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

### Construct sequence table #### 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

### Remove chimeras ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

### Track reads through the pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

### Assign taxonomy ####

## HERE is where you can add the CyanoSeq database to better classify the cyano sequences ##
# put the fasta file with the silva classification in the path folder

taxa <- assignTaxonomy(seqtab.nochim, "/Users/briannepalmer/Downloads/CyanoSeq_1.1.1_dada2.fastq.gz", multithread=TRUE); beep(sound = 3)
#taxa <- addSpecies(taxa, "/Users/briannepalmer/Downloads/silva_species_assignment_v138.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)



### Handoff to phyloseq ####

library(phyloseq)
library(Biostrings)
library(ggplot2)
library(tidyverse)
theme_set(theme_bw())

# We can construct a simple sample data.frame from the information encoded in the filenames. Usually this step would instead involve reading the sample data in from a file.

# change this to match your samples 

samples.out <- rownames(seqtab.nochim)
sampleID <- c("R18", "R19", "R20", "R22")
Treatment <- c("BSC_600_2cm", "BSC_600_1cm", "BSC_600_8cm", "BSC_Cont_5cm")
Type <- c("BSC", "BSC", "BSC", "BSC")
Temp <- c("600", "600", "600", "Cont")
Depth<- c("2cm", "1cm", "8cm", "5cm")

metadata <- cbind(sampleID, Type, Temp, Depth, Treatment)

samdf <- as.data.frame(metadata)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

taxa.dna <- data.frame(ps@tax_table)
cyano.sequences <- taxa.dna %>% dplyr::filter(Phylum == "Cyanobacteriota")
cyano.sequences$Sequence <- rownames(cyano.sequences)

fwrite(cyano.sequences, "/Users/briannepalmer/Desktop/cyano_sequences.txt", sep = ",")


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

## make OTU table for other analyses ####
#taxa_names(ps) <- paste0("SV", seq(ntaxa(ps)))
# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(ps), "matrix")
# transpose if necessary
#if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

# make tax table based on relative abundance ####
TGroup <- tax_glom(ps, taxrank = "Genus")

# write OTU and TAX tables ####
fwrite(OTUdf, "/Users/briannepalmer/Desktop/otu_cy_fire_cyanoSeq_label.txt", sep = ",")
write.csv(ps@tax_table, "/Users/briannepalmer/Desktop/tax_cy_fire_cyanoSEQ.csv")

otu <- data.table::fread(file = "/Users/briannepalmer/Desktop/Chapter 4 Trees/otu_cy_fire_cyanoSeq.txt")
taxa <- data.table::fread(file ="/Users/briannepalmer/Desktop/Chapter 4 Trees/tax_cy_fire_cyanoSeq.csv", header = TRUE)
metadata <- as.matrix(metadata)

otu.t <- t(otu)
colnames(otu.t) <- metadata[,1]
V1 <- rownames(otu.t)
otu.t <- cbind(otu.t, V1)
otu.t <- unlist(otu.t)
otu.m <- merge(otu.t, taxa, by = "V1")

#look for sequences that are cyanobacteria and put them in a new text file 
#cyano <- otu.m %>% filter(Phylum == "Cyanobacteriota")

#fwrite(cyano, "/Users/briannepalmer/Desktop/all_cyano02.txt", sep = ",")

otu.m[c(2:5)] <- sapply(otu.m[c(2:5)],as.numeric)

data.t <- otu.m %>% pivot_longer(cols = R18:R22, names_to = "Sample", values_to = "count")
all.data <- cbind(metadata, data.t)


plot_richness(ps, x="Treatment", measures=c("Shannon", "Simpson"))
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Treatment", title="Bray NMDS")

library(ecodist)
library(vegan)
#make bray-curtis matrix
otu.bray <- vegdist(ps.prop@otu_table, method = "bray")# dissimilarity matrix using bray-curtis distance indices
#Calculate PCoA values 
pcoaVS <- pco(otu.bray, negvals = "zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero

pcoa1 = pcoaVS$vectors[,1] # adds the points of the NMDS 1 dimmension
pcoa2 = pcoaVS$vectors[,2] # adds the points of the NMDS 2 dimmension

MDS = data.frame(pcoa1 = pcoa1, pcoa2 = pcoa2)
MDS = cbind(MDS, metadata)

smallpal2 <- palette(c("#0073C2FF", "#EFC000FF" , "#f08080", "#CD534CFF"))
smallpal2 <- colorRampPalette(smallpal2)(5)
ggplot(MDS, aes(x = pcoa1, y = pcoa2, color = Treatment)) + geom_point(size =3) + 
  theme_bw() + stat_ellipse(aes(color = Treatment)) + scale_color_manual(values = smallpal2)

metadata <- data.frame(metadata)

otu.table <- data.frame(ps.prop@otu_table)
adonis2(otu.table ~ Treatment, metadata) 


phylum <-  all.data %>% group_by(Treatment, Phylum) %>% summarize(n = sum(count))%>% mutate(rel = n / sum(n))
phylum <- phylum %>% dplyr::select(Treatment, Phylum, rel)
phylum <- phylum %>% pivot_wider(names_from = Phylum, values_from = rel)
phylum[c(2:7)] <- sapply(phylum[c(2:7)],as.numeric)
phylum[is.na(phylum)] <- 0

phylum <- as.data.frame(c(metadata, phylum))
phylum.ggplot <- pivot_longer(phylum, cols = Actinobacteria:NA., names_to = "phylum", values_to = "count")

ggplot(phylum.ggplot, aes(fill = phylum, y = count, x = Treatment)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance") 

ggplot(phylum.ggplot %>% filter(phylum != "NA."), aes(fill = phylum, y = count, x = Treatment)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme_bw()

phylum.count <- phylum.ggplot %>%  group_by(phylum) %>% summarize(sum = sum(count))

sum(phylum.count$sum)
2.6043038917/4


t.test(count ~ sampleID, all.data %>% filter(phylum == "Cyanobacteriota") %>% filter(sampleID == c("R18", "R19")))

no.na <- all.data %>% filter(Phylum != "NA.")
no.na <- no.na %>% pivot_wider(names_from = V1, values_from = count)
no.na <- no.na %>% 
  mutate_at(c(13:366), as.numeric) 
no.na[c(13:366)][is.na(no.na[c(13:366)])]<- 0
no.na$sum <- rowSums(no.na[c(13:366)])
no.na <- no.na %>% filter(sum >0)
otu.no.na <- no.na[c(13:366)]
otu.no.na <- as.numeric(otu.no.na)
otu.no.na[is.na(otu.no.na)]<- 0
md <- no.na[c(1:12)]


all.dat <- all.data %>% pivot_wider(names_from = V1, values_from = count)
all.dat <- all.dat %>% 
  mutate_at(c(13:874), as.numeric) 
all.dat[c(13:874)][is.na(all.dat[c(13:874)])]<- 0
all.dat$sum <- rowSums(all.dat[c(13:874)])
all.dat <- all.dat %>% filter(sum >0)
otu.all.dat <- all.dat[c(13:874)]
otu.all.dat <- as.numeric(otu.all.dat)
otu.all.dat[is.na(otu.all.dat)]<- 0
md <- all.dat[c(1:12)]

adonis2(otu.all.dat ~ Treatment, md) 

proteo <-  all.data %>% filter(Phylum == "Proteobacteria") %>% group_by(Treatment, Class) %>% summarize(n = sum(count))%>% mutate(rel = n / sum(n)) 
proteo <- proteo %>% dplyr::select(Treatment, Class, rel)
proteo <- proteo %>% pivot_wider(names_from = Class, values_from = rel)
proteo[c(2:5)] <- sapply(proteo[c(2:5)],as.numeric)
proteo[is.na(proteo)] <- 0

proteo <- as.data.frame(c(metadata, proteo))
proteo.ggplot <- pivot_longer(proteo, cols = Alphaproteobacteria:NA., names_to = "Class", values_to = "count")

ggplot(proteo.ggplot, aes(fill = Class, y = count, x = Treatment)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance") 

proteo <-  all.data %>% filter(Phylum == "Proteobacteria") %>% group_by(Treatment, Order) %>% summarize(n = sum(count))%>% mutate(rel = n / sum(n)) 
proteo <- proteo %>% dplyr::select(Treatment, Order, rel)
proteo <- proteo %>% pivot_wider(names_from = Order, values_from = rel)
proteo[c(2:7)] <- sapply(proteo[c(2:7)],as.numeric)
proteo[is.na(proteo)] <- 0

proteo <- as.data.frame(c(metadata, proteo))
proteo.ggplot <- pivot_longer(proteo, cols = Burkholderiales:NA., names_to = "Order", values_to = "count")

ggplot(proteo.ggplot, aes(fill = Order, y = count, x = Treatment)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance") 



# there are a lot of NAs that could be identified using the silva database, it could be worth comparing the two databases or combining them. 


############################

taxa <- assignTaxonomy(seqtab.nochim, "/Users/briannepalmer/Downloads/silva_nr_v138_train_set.fa.gz", multithread=TRUE); beep(sound = 3)
#taxa <- addSpecies(taxa, "/Users/briannepalmer/Downloads/silva_species_assignment_v138.fa.gz")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# We can construct a simple sample data.frame from the information encoded in the filenames. Usually this step would instead involve reading the sample data in from a file.

# change this to match your samples 

samples.out <- rownames(seqtab.nochim)
sampleID <- c("R18", "R19", "R20", "R22")
Treatment <- c("BSC_600_2cm", "BSC_600_1cm", "BSC_600_8cm", "BSC_Cont_5cm")
Type <- c("BSC", "BSC", "BSC", "BSC")
Temp <- c("600", "600", "600", "Cont")
Depth<- c("2cm", "1cm", "8cm", "5cm")

metadata <- cbind(sampleID, Type, Temp, Depth, Treatment)

samdf <- as.data.frame(metadata)
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

## make OTU table for other analyses ####
#taxa_names(ps) <- paste0("SV", seq(ntaxa(ps)))
# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(ps), "matrix")
# transpose if necessary
#if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

# make tax table based on relative abundance ####
TGroup <- tax_glom(ps, taxrank = "Genus")

# write OTU and TAX tables ####
fwrite(OTUdf, "/Users/briannepalmer/Desktop/otu_cy_fire_SILVA.txt", sep = ",")
write.csv(ps@tax_table, "/Users/briannepalmer/Desktop/tax_cy_fire_SILVA.csv")

otu <- data.table::fread(file = "/Users/briannepalmer/Desktop/Chapter 4 Trees/otu_cy_fire_SILVA.txt")
taxa <- data.table::fread(file ="/Users/briannepalmer/Desktop/Chapter 4 Trees/tax_cy_fire_SILVA.csv", header = TRUE)

otu.t <- t(otu)
colnames(otu.t) <- metadata[,1]
V1 <- rownames(otu.t)
otu.t <- cbind(otu.t, V1)
otu.t <- unlist(otu.t)
otu.m <- merge(otu.t, taxa, by = "V1")

# look for sequences that are cyanobacteria and put them in a new text file 
#cyano <- otu.m %>% filter(Phylum == "Cyanobacteriota")
#fwrite(cyano, "Desktop/all_cyano.txt", sep = ",")

otu.m[c(2:5)] <- sapply(otu.m[c(2:5)],as.numeric)

data.t <- otu.m %>% pivot_longer(cols = R18:R22, names_to = "Sample", values_to = "count")
all.data <- cbind(metadata, data.t)


plot_richness(ps, x="Treatment", measures=c("Shannon", "Simpson"))
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Treatment", title="Bray NMDS")

library(ecodist)
library(vegan)
#make bray-curtis matrix
otu.bray <- vegdist(ps.prop@otu_table, method = "bray")# dissimilarity matrix using bray-curtis distance indices
#Calculate PCoA values 
pcoaVS <- pco(otu.bray, negvals = "zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero

pcoa1 = pcoaVS$vectors[,1] # adds the points of the NMDS 1 dimmension
pcoa2 = pcoaVS$vectors[,2] # adds the points of the NMDS 2 dimmension

MDS = data.frame(pcoa1 = pcoa1, pcoa2 = pcoa2)
MDS = cbind(MDS, metadata)

smallpal2 <- palette(c("#0073C2FF", "#EFC000FF" , "#f08080", "#CD534CFF"))
smallpal2 <- colorRampPalette(smallpal2)(5)
ggplot(MDS, aes(x = pcoa1, y = pcoa2, color = Treatment)) + geom_point(size =3) + 
  theme_bw() + stat_ellipse(aes(color = Treatment)) + scale_color_manual(values = smallpal2)

metadata <- data.frame(metadata)

otu.table <- data.frame(ps.prop@otu_table)
adonis2(otu.table ~ Treatment, metadata) # p = 0.001

library(pairwiseAdonis)
set.seed(12345)
pairwise.adonis(otu.table, metadata$Treatment, p.adjust.m = "BH") 

phylum <-  all.data %>% group_by(Treatment, Phylum) %>% summarize(n = sum(count))%>% mutate(rel = n / sum(n))
phylum <- phylum %>% dplyr::select(Treatment, Phylum, rel)
phylum <- phylum %>% pivot_wider(names_from = Phylum, values_from = rel)
phylum[c(2:7)] <- sapply(phylum[c(2:7)],as.numeric)
phylum[is.na(phylum)] <- 0

phylum <- as.data.frame(c(metadata, phylum))
phylum.ggplot <- pivot_longer(phylum, cols = Acidobacteriota:NA., names_to = "phylum", values_to = "count")

ggplot(phylum.ggplot, aes(fill = phylum, y = count, x = Treatment)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance") 

ggplot(phylum.ggplot %>% filter(phylum != "NA."), aes(fill = phylum, y = count, x = Treatment)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance")+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme_bw()

phylum.count <- phylum.ggplot %>%  group_by(phylum) %>% summarize(sum = sum(count))

sum(phylum.count$sum)

3.3065235782/4

no.na <- all.data %>% filter(Phylum != "NA.")
no.na <- no.na %>% pivot_wider(names_from = V1, values_from = count)
no.na <- no.na %>% 
  mutate_at(c(13:286), as.numeric) 
no.na[c(13:286)][is.na(no.na[c(13:286)])]<- 0
no.na$sum <- rowSums(no.na[c(13:286)])
no.na <- no.na %>% filter(sum >0)
otu.no.na <- no.na[c(13:286)]
otu.no.na <- as.numeric(otu.no.na)
otu.no.na[is.na(otu.no.na)]<- 0
md <- no.na[c(1:12)]


all.dat <- all.data %>% pivot_wider(names_from = V1, values_from = count)
all.dat <- all.dat %>% 
  mutate_at(c(13:913), as.numeric) 
all.dat[c(13:913)][is.na(all.dat[c(13:913)])]<- 0
all.dat$sum <- rowSums(all.dat[c(13:913)])
all.dat <- all.dat %>% filter(sum >0)
otu.all.dat <- all.dat[c(13:913)]
otu.all.dat <- as.numeric(otu.all.dat)
otu.all.dat[is.na(otu.all.dat)]<- 0
md <- all.dat[c(1:12)]

adonis2(otu.all.dat ~ Treatment, md) # 


proteo <-  all.data %>% filter(Phylum == "Proteobacteria") %>% group_by(Treatment, Class) %>% summarize(n = sum(count))%>% mutate(rel = n / sum(n)) 
proteo <- proteo %>% dplyr::select(Treatment, Class, rel)
proteo <- proteo %>% pivot_wider(names_from = Class, values_from = rel)
proteo[c(2:4)] <- sapply(proteo[c(2:4)],as.numeric)
proteo[is.na(proteo)] <- 0

proteo <- as.data.frame(c(metadata, proteo))
proteo.ggplot <- pivot_longer(proteo, cols = Alphaproteobacteria:NA., names_to = "Class", values_to = "count")

ggplot(proteo.ggplot, aes(fill = Class, y = count, x = Treatment)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance") 

proteo <-  all.data %>% filter(Phylum == "Proteobacteria") %>% group_by(Treatment, Order) %>% summarize(n = sum(count))%>% mutate(rel = n / sum(n)) 
proteo <- proteo %>% dplyr::select(Treatment, Order, rel)
proteo <- proteo %>% pivot_wider(names_from = Order, values_from = rel)
proteo[c(2:7)] <- sapply(proteo[c(2:7)],as.numeric)
proteo[is.na(proteo)] <- 0

proteo <- as.data.frame(c(metadata, proteo))
proteo.ggplot <- pivot_longer(proteo, cols = Acetobacterales:NA., names_to = "Order", values_to = "count")

ggplot(proteo.ggplot, aes(fill = Order, y = count, x = Treatment)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance") 


