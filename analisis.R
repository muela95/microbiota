library(dplyr)
library(dada2)
library(readxl)
library(ggplot2)
library(msa)
library(phyloseq)
library(haplotypes)
library(phangorn)
library(microbiome)
library(ggpubr)
library(vegan)
library(tidyr)
library(metafolio)
library(ggforce)
library(ggbeeswarm)
library(tibble)
library(Tjazi)
library(stringr)
library(rstatix)
library(kableExtra)

set.seed(123)
setwd("D:/Laboratorio/Microbiota yo/2024")
metadata <- read_excel("metadata.xlsx")
metadata$Group <- factor(metadata$Group, levels = c("Sed", "Run", "RunTime", 
                                                    "RunVel"))

conducta <- read_excel("D:/Laboratorio/Microbiota y conducta/Excel Exploratorio Microbiota Nivel1 Reino.xlsx", 
                       sheet = "Hoja1")
metadata$CI <- conducta$IDPromedio
paletaEli <- c("#a0a0a4", "#0f99b2", "#05be78", "#056943")


path <- ("D:/Laboratorio/Microbiota yo/2023/ranomaly/raw")
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,210),
#                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#                    compress=TRUE, multithread=F, trimLeft = c(17,21)) 

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250),
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=F, trimLeft = c(37)) 

head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
#Este ultimo comando te da la frecuencia de secuencias de esa longitud
#Si quisieses quedarte con un rango de longitudes (rollo de 439 a 465)
#Puedes hacer un seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(439,465)]
#y ojo que son secuencias en número de bases, no de pares de bases

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")
rownames(track) <- sample.names
head(track)
track <- as.data.frame(track)
track$porcentaje <- track$nonchim/track$input
write.table(track, file="lecturas.tsv", quote = F, sep="\t", col.names=NA)

taxa <- assignTaxonomy(seqtab.nochim, "D:/Laboratorio/Microbiota yo/2023/ranomaly/silva_nr99_v138.1_train_set.fa.gz", multithread=T)
taxa <- addSpecies(taxa, "D:/Laboratorio/Microbiota yo/2023/ranomaly/silva_species_assignment_v138.1.fa.gz")
#creo que también se podría hacer del tirón con solo un comando con el otro 
#archivo que es todo ya con especies

samples.out <- rownames(seqtab.nochim)
metadata <- as.data.frame(metadata[2])
rownames(metadata) <- samples.out

seqs <- getSequences(seqtab)
names(seqs) <- seqs 
mult <- msa(seqs, method="ClustalW", type="dna", order="input")

phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", 
                    control = pml.control(trace = 0))

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               phyloseq::sample_data(metadata), 
               phyloseq::tax_table(taxa),
               phyloseq::phy_tree(fitGTR$tree))


saveRDS(ps, file = "phyloseq_object.rds")










