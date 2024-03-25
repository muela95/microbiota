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





#########################################################################
# New code starts here (besides changes in filtering, everything 
# above is the same)

# First step is to adjust our data to look like Thomaz's; samples on columns,
# bacteria on rows. Since our Genus level has quite a few NA values, we will
# work with Genus and Family level too.

# Everything needs to be a dataframe in order for next part of code to work
taxonomia <- as.data.frame(tax_table(ps))
generos <- as.data.frame(otu_table(ps))
familias <- as.data.frame(otu_table(ps))

# taxonomia, generos and familias all three are in the same order (rows to
# columns, respectively) so no need of search and replace sequences.
colnames(familias) <- taxonomia$Family
colnames(generos) <- taxonomia$Genus


# Transposing count tables so we have samples as columns:
familias <- t(familias)
generos <- t(generos)
# Apparently t() turns dataframes to matrices, so again, we need to turn them to
# dataframes
 
familias <- as.data.frame(familias)
generos <- as.data.frame(generos)

# Repeating rownames now have a suffix ".1", ".2", ".3"... which is not ideal,
# so, function to detect anything repeating before a dot in rownames, and add 
# all cells of that row. 

summarize_by_prefix <- function(df) {
  df <- df %>%
    mutate(prefix = sub("\\..*", "", rownames(df)))
  df_summary <- df %>%
    group_by(prefix) %>%
    summarize_all(sum)
  return(df_summary)
}


countFamilias <- summarize_by_prefix(familias)
countGeneros <- summarize_by_prefix(generos)
# This column of metadata needed to be the same
metadata$sample <- colnames(countGeneros)[2:20]

# countFamilias and countGeneros need to be dataframes, and take values of a 
# column (prefix) as rownames
countFamilias <- as.data.frame(countFamilias)
rownames(countFamilias) <- countFamilias$prefix
countFamilias$prefix <- NULL # We won't need that column anymore
countGeneros <- as.data.frame(countGeneros)
rownames(countGeneros) <- countGeneros$prefix
countGeneros$prefix <- NULL

# Making sure every value is numeric
countGeneros   <- apply(countGeneros,c(1,2),function(x) 
  as.numeric(as.character(x)))
countFamilias   <- apply(countFamilias,c(1,2),function(x) 
  as.numeric(as.character(x)))

# There's NAs, so we'll remove that row
countFamilias <- countFamilias[-which(rownames(countFamilias) == "NA"),]
countGeneros <- countGeneros[-which(rownames(countGeneros) == "NA"),]


# Making an untouched, unfiltered copy now, overwriting previous dataframe that
# we'll no longer need. Also, exporting .csv
familias <- as.data.frame(countFamilias)
generos <- as.data.frame(countGeneros)
write.csv(familias, "D:/Laboratorio/Elisa/Microbiota2024/familias.csv",
                     row.names=T)
write.csv(generos, "D:/Laboratorio/Elisa/Microbiota2024/generos.csv",
                     row.names=T)
write.csv(taxonomia, "D:/Laboratorio/Elisa/Microbiota2024/taxonomia.csv",
          row.names=T)


# We need to remove bacteria with a prevalence fo less than 10% (each bacteria 
# has to be in at least 10% of all samples). First we calculate, for each 
# bacteria, in how many samples they are not present.
n_zeroesG <- rowSums(countGeneros == 0)
n_zeroesF <- rowSums(countFamilias == 0)

# Now we keep only rows with a number of zeroes less or equal to 90% of the
# total number of samples
countGeneros <- countGeneros[n_zeroesG <= round(ncol(countGeneros)*0.9),]
countFamilias <- countFamilias[n_zeroesF <= round(ncol(countFamilias)*0.9),]

# To check how many we removed:
cat("We removed" ,length(n_zeroesF)-dim(countFamilias)[1], "from Familias")
cat("We removed" ,length(n_zeroesG)-dim(countGeneros)[1], "from Generos")

# Centered Log-Ratio transform to the dataframes; log(value/geometric mean).
# Zeroes are replaced with 65% of the lowest value
generosCLR <- clr_c(countGeneros)
familiasCLR <- clr_c(countFamilias)


############################################ Stacked barplots
bargenus <- generos
barfamilies <- familias


#Make into relative abundance
bargenus   <- apply(bargenus, 2, function(i) i/sum(i)) 
barfamilies   <- apply(barfamilies, 2, function(i) i/sum(i)) 

#Define a cutoff for rare taxa in several steps:
#first, determine the max % abundance every feature ever shows up at 
maxabundancesg <- apply(bargenus, 1, max)
maxabundancesf <- apply(barfamilies, 1, max)
#Meanwhile, transpose the count table for future wrangling.
bargenus      <- data.frame(t(bargenus))
barfamilies      <- data.frame(t(barfamilies))

#For every sample, sum up all rare taxa ( < 1% at their highest in this case)
bargenus$`Rare Taxa` <- rowSums(bargenus[,maxabundancesg < 0.01], na.rm = TRUE)
barfamilies$`Rare Taxa` <- rowSums(barfamilies[,maxabundancesf < 0.01], 
                                   na.rm = TRUE)

#Remove the individual rare taxa now that they're summed up
bargenus = bargenus[,c(maxabundancesg > 0.01, T) ] #`T` to include Rare Taxa  
barfamilies = barfamilies[,c(maxabundancesf > 0.01, T) ] 
#Prepare the data for ggplot by adding in metadata here
bargenus$Group       = metadata$Group
bargenus$ID          = metadata$sample
barfamilies$Group       = metadata$Group
barfamilies$ID          = metadata$sample

#Wrangle the data to long format for easy plotting
barlongg = bargenus %>% 
  pivot_longer(!c(ID, Group), names_to = c("Microbe"), values_to = "value") %>%
  mutate(Microbe = str_replace(Microbe, ".*_or_", ""))
barlongf = barfamilies %>% 
  pivot_longer(!c(ID, Group), names_to = c("Microbe"), values_to = "value") %>%
  mutate(Microbe = str_replace(Microbe, ".*_or_", ""))

#Change the colour for the rare taxa to gray to make them stand out
colsg = metafolio::gg_color_hue(length(unique(barlongg$Microbe)))
colsg[unique(barlongg$Microbe)=="Rare Taxa"]="dark gray"
colsf = metafolio::gg_color_hue(length(unique(barlongf$Microbe)))
colsf[unique(barlongf$Microbe)=="Rare Taxa"]="dark gray"

#Create the stacked barplots using ggplot2
barlongg %>%
  ggplot(aes(x = ID, y = value, fill = Microbe)) + 
  geom_bar(stat = "identity", col = "black", size = .2, width = 1) + 
  facet_row(~Group, scales = "free_x") +
  #Adjust layout and appearance
  scale_fill_manual(values = colsg, labels = unique(sub(".*ales_", "", 
                                                       barlongg$Microbe))) + 
  ggtitle("Genus")+
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 1, keyheight = 1, title = "Legend")) + 
  theme_bw() + xlab("") +  ylab("Proportion") + 
  theme(text = element_text(size = 14), axis.text.x = element_blank())

barlongf %>%
  ggplot(aes(x = ID, y = value, fill = Microbe)) + 
  geom_bar(stat = "identity", col = "black", size = .2, width = 1) + 
  facet_row(~Group, scales = "free_x") +
  #Adjust layout and appearance
  scale_fill_manual(values = colsf, labels = unique(sub(".*ales_", "", 
                                                       barlongf$Microbe))) + 
  ggtitle("Families")+
  scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(ncol = 1, keyheight = 1, title = "Legend")) + 
  theme_bw() + xlab("") +  ylab("Proportion") + 
  theme(text = element_text(size = 14), axis.text.x = element_blank())


# Exporto manualmente las dos a 1000x1000 en .svg



################################# Alpha diversity
alpha_genus = get_asymptotic_alpha(species = generos, verbose = F)
alpha_families = get_asymptotic_alpha(species = familias, verbose = F)

alpha_genus$Legend = metadata$Group
alpha_genus$ID    = metadata$sample
alpha_families$Legend = metadata$Group
alpha_families$ID    = metadata$sample



alpha_genus %>%
  pivot_longer(!c(Legend, ID)) %>%
  ggplot(aes(x     = Legend,
             y     = value, 
             fill  = Legend)) + 
  geom_boxplot(alpha = 1/2, coef = 100) + 
  geom_jitter(position = position_jitterdodge()) + 
  facet_wrap(~name, scales = "free") + theme_bw()  +
  scale_fill_manual(values = paletaEli)  +
  stat_compare_means(method="anova", label.y.npc = "bottom") +
  ylab("") + xlab("") + ggtitle("Genus")

alpha_families %>%
  pivot_longer(!c(Legend, ID)) %>%
  ggplot(aes(x     = Legend,
             y     = value, 
             fill  = Legend)) + 
  geom_boxplot(alpha = 1/2, coef = 100) + 
  geom_jitter(position = position_jitterdodge()) + 
  facet_wrap(~name, scales = "free") + theme_bw()  +
  scale_fill_manual(values = paletaEli)  +
  stat_compare_means(method="anova", label.y.npc = "bottom") +
  ylab("") + xlab("") + ggtitle("Family")


anova_test(alpha_genus, Chao1 ~ Legend)                # p=0.03
tukey_hsd(alpha_genus, Chao1 ~ Legend)                 # Run-RunVel p=0.027
anova_test(alpha_genus, `Shannon Entropy` ~ Legend)    # p=0.051
anova_test(alpha_genus, `Simpson Index` ~ Legend)      # p=0.071

anova_test(alpha_families, Chao1 ~ Legend)             # n.s.
anova_test(alpha_families, `Shannon Entropy` ~ Legend) # n.s.
anova_test(alpha_families, `Simpson Index` ~ Legend)   # n.s.




############################################## Beta diversity

pcaGenus <- prcomp(t(generosCLR))
pcaFamilies <- prcomp(t(familiasCLR))

pc1g <- round(pcaGenus$sdev[1]^2/sum(pcaGenus$sdev^2),4) * 100
pc2g <- round(pcaGenus$sdev[2]^2/sum(pcaGenus$sdev^2),4) * 100
pc1f <- round(pcaFamilies$sdev[1]^2/sum(pcaFamilies$sdev^2),4) * 100
pc2f <- round(pcaFamilies$sdev[2]^2/sum(pcaFamilies$sdev^2),4) * 100

pcaF <- data.frame(PC1 = pcaFamilies$x[,1],
                   PC2 = pcaFamilies$x[,2])
pcaF$ID <- metadata$Animal
pcaF$Legend <- metadata$Group

pcaG <- data.frame(PC1 = pcaGenus$x[,1],
                   PC2 = pcaGenus$x[,2])
pcaG$ID <- metadata$Animal
pcaG$Legend <- metadata$Group


ggplot(pcaF, aes(x=PC1, y=PC2, fill=Legend, color=Legend, group=Legend))+
  stat_ellipse(geom="polygon", alpha=0.25) + geom_point(size=3)+
  scale_fill_manual(values=paletaEli) + scale_color_manual(values=paletaEli)+
  ggtitle("Family") + xlab(paste("PC1: ", pc1f, "%", sep="")) +
  ylab(paste("PC2: ", pc2f, "%", sep="")) + theme_bw()

ggplot(pcaG, aes(x=PC1, y=PC2, fill=Legend, color=Legend, group=Legend))+
  stat_ellipse(geom="polygon", alpha=0.25) + geom_point(size=3)+
  scale_fill_manual(values=paletaEli) + scale_color_manual(values=paletaEli)+
  guides(fill=guide_legend(override.aes = list(shape=c(21)))) + 
  ggtitle("Genus") + xlab(paste("PC1: ", pc1g, "%", sep="")) +
  ylab(paste("PC2: ", pc2g, "%", sep="")) + theme_bw()


options(knitr.kable.NA = "")
dist_fam = dist(t(familiasCLR), method="euclidean")
dist_gen = dist(t(generosCLR), method="euclidean")
beta_disp_fam = betadisper(dist_fam, group = metadata$Group)
beta_disp_gen = betadisper(dist_gen, group = metadata$Group)

kable(anova(beta_disp_fam), digits = 4)
kable(anova(beta_disp_gen), digits = 4)

permanova_fam <- adonis2(dist_fam ~ Group, data = metadata, 
                         method = "euclidean", permutations = 1000)
kable(permanova_fam, digits=4)

permanova_gen <- adonis2(dist_gen ~ Group, data = metadata,
                         method = "euclidean", permutations = 1000)
kable(permanova_gen, digits=4)


############################################### Differential abundance


genus.glm = fw_glm(x             = generosCLR,
                   f             = ~ Group, 
                   metadata      = metadata, 
                   adjust.method = "BH", format = "brief")
family.glm = fw_glm(x             = familiasCLR,
                   f             = ~ Group, 
                   metadata      = metadata, 
                   adjust.method = "BH", format = "brief")
glimpse(genus.glm)

hist(genus.glm$`GroupRun Pr(>|t|)`, xlim = c(0, 1), breaks = 20)
# OJO, bimodal??

# fw_glm() seem to compare every group using Sed as reference

significant <- c()

todoConTodo <- function(df){
  tdf <- as.data.frame(t(df))
  tdf$Group <- metadata$Group
  for(i in colnames(tdf)[1:(ncol(tdf) - 1)]) {
    model <- aov(as.formula(paste(i, " ~ Group")), data = tdf)
    summarymodel <- summary(model)
    if(summarymodel[[1]]$`Pr(>F)`[1] <= 0.05){
      significant <- c(significant, i)
      print(i)
      tukey <- TukeyHSD(model)
      print(tukey)
    }
  }
  return(significant)
}


significativosFam <- todoConTodo(familiasCLR)
significativosGen <- todoConTodo(generosCLR)

# Every ANOVA with a p <= 0.05 in families
familiasCLR %>%
  t() %>%
  as.data.frame()  %>%
  add_column(Group = metadata$Group) %>%
  pivot_longer(!c("Group"))  %>%
  filter(name %in% significativosFam) %>%
  ggplot(aes(x     = Group, 
             y     = value, 
             fill  = Group, 
             group = Group)) + 
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_beeswarm(size = 3, cex = 3) + 
  facet_wrap(~name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = paletaEli) + 
  stat_compare_means(method="anova") +
  ggtitle("Families") +
  ylab("") + xlab("") + theme_bw() + theme(text = element_text(size = 12))

# Every ANOVA with a p <= 0.05 in genus
generosCLR %>%
  t() %>%
  as.data.frame()  %>%
  add_column(Group = metadata$Group) %>%
  pivot_longer(!c("Group"))  %>%
  filter(name %in% significativosGen) %>%
  ggplot(aes(x     = Group, 
             y     = value, 
             fill  = Group, 
             group = Group)) + 
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_beeswarm(size = 3, cex = 3) + 
  facet_wrap(~name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = paletaEli) + 
  stat_compare_means(method="anova") +
  ggtitle("Genus") +
  ylab("") + xlab("") + theme_bw() + theme(text = element_text(size = 12))
  
  
  
  
# Every differentially expressed bacteria at family level correlated with 
# behavior
familiasCLR %>%
  t() %>%
  as.data.frame()  %>%
  add_column(Group = metadata$Group, .before = 1) %>%
  add_column(CI = metadata$CI, .after = "Group") %>%
  pivot_longer(!c("Group", "CI"))  %>% 
  filter(name %in% significativosFam)  %>% 
  ggscatter("value", "CI", 
            color = "Group",fill="lightgrey", palette=paletaEli) + 
  stat_cor(method = "pearson")+
  facet_wrap(~name, scales="free_x") + geom_hline(yintercept = 0.2, linetype=
                                                    "dashed") +
  geom_hline(yintercept = 0) + geom_smooth(method="lm", color="black") +
  ggtitle("Families")

  
  

# Every differentially expressed bacteria at genus level correlated with 
# behavior
generosCLR %>%
  t() %>%
  as.data.frame()  %>%
  add_column(Group = metadata$Group, .before = 1) %>%
  add_column(CI = metadata$CI, .after = "Group") %>%
  pivot_longer(!c("Group", "CI"))  %>% 
  filter(name %in% significativosGen)  %>% 
  ggscatter("value", "CI", 
            color = "Group",fill="lightgrey", palette=paletaEli) + 
  stat_cor(method = "pearson")+
  facet_wrap(~name, scales="free_x") + geom_hline(yintercept = 0.2, linetype=
                                                    "dashed") +
  geom_hline(yintercept = 0) + geom_smooth(method="lm", color="black")






# Correlation matrices
corFam <- familiasCLR %>%
  t() %>%
  as.data.frame()  %>%
  add_column(CI = metadata$CI, .before = 1) 

corFam <- corFam[,c(significativosFam, "CI")]
View(corFam)
rhoFam <- cor(corFam)
pcorFam <- cor_pmat(corFam)
ggcorrplot(rhoFam, p.mat = pcorFam, type = "lower", insig = "blank", lab=T, 
           colors = c("#a23203", "white", "#0373a2"))



corGen <- generosCLR %>%
  t() %>%
  as.data.frame()  %>%
  add_column(CI = metadata$CI, .before = 1) 

corGen <- corGen[,c(significativosGen, "CI")]
rhoGen <- cor(corGen)
pcorGen <- cor_pmat(corGen)
ggcorrplot(rhoGen, p.mat = pcorGen, type = "lower", insig = "blank", lab=T, 
           colors = c("#a23203", "white", "#0373a2"))





corSigGen <- significativosGen[c(2,3,5,9,11)]


generosCLR %>%
  t() %>%
  as.data.frame()  %>%
  add_column(Group = metadata$Group, .before = 1) %>%
  add_column(CI = metadata$CI, .after = "Group") %>%
  pivot_longer(!c("Group", "CI"))  %>% 
  filter(name %in% corSigGen)  %>% 
  ggscatter("value", "CI", 
            color = "Group",fill="lightgrey", palette=paletaEli) + 
  stat_cor(method = "pearson")+
  facet_wrap(~name, scales="free_x", ncol = 5) + geom_hline(yintercept = 0.2, 
                                                            color="red") +
  geom_hline(yintercept = 0) + geom_smooth(method="lm", color="black")




psSedRun <- phyloseq::subset_samples(
  enterotypes_arumugam,
  Group %in% c("Sed", "Run")
)

psSedRV <- phyloseq::subset_samples(
  enterotypes_arumugam,
  Group %in% c("Sed", "RunVel")
)

psSedRT <- phyloseq::subset_samples(
  enterotypes_arumugam,
  Group %in% c("Sed", "RunTime")
)

psRunRT <- phyloseq::subset_samples(
  enterotypes_arumugam,
  Group %in% c("Run", "RunTime")
)

psRunRV <- phyloseq::subset_samples(
  enterotypes_arumugam,
  Group %in% c("Run", "RunVel")
)

psRTRV <- phyloseq::subset_samples(
  enterotypes_arumugam,
  Group %in% c("RunTime", "RunVel")
)

run_aldex(psSedRun, group="Group")















