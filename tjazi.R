library(phyloseq)
library(Tjazi)
library(readxl)
library(ggcorrplot)
library(rstatix)
library(tidyr)
library(stringr)
library(ggforce)
library(ggpubr)
library(vegan)
library(kableExtra)
library(tibble)
library(MicrobiotaProcess)
library(UpSetR)
library(patchwork)

set.seed(1)

paletaEli <- c("#a0a0a4", "#0f99b2", "#05be78", "#056943")

# Set working directory and load phyloseq object (check analisis.R to see 
# filtering and crop parameters, updated after Thomaz's suggestions).
setwd("D:/Laboratorio/Elisa/Microbiota2024")
ps <- readRDS("phyloseq_object.rds")
metadata <- read_excel("metadata.xlsx")
metadata$Group <- factor(metadata$Group, levels = c("Sed", "Run", "RunTime", 
                                                    "RunVel"))


# We can merge our object at any taxonomic rank with "tax_glom()", since we may
# have a low number of identified Genus, we'll be working with family level too,
# just in case.
psG <- tax_glom(ps, taxrank = "Genus")
psF <- tax_glom(ps, taxrank = "Family")



taxonomiaG <- as.data.frame(tax_table(psG))
taxonomiaF <- as.data.frame(tax_table(psF))
generos <- as.data.frame(otu_table(psG))
familias <- as.data.frame(otu_table(psF))

# Tax tables and "otu" tables have the same order within the same phyloseq
# object, so we can just replace colnames instead of iterating, searching and
# replacing individual values.
# This can be checked with:
# colnames(familias) == rownames(taxonomiaF)
# and
# colnames(generos) == rownames(taxonomiaG)
colnames(familias) <- taxonomiaF$Family
colnames(generos) <- taxonomiaG$Genus


# Transposing count tables so we have samples as columns:
familias <- t(familias)
generos <- t(generos)
# "t()" turns dataframes to matrices, so again, we need to turn them to
# dataframes

familias <- as.data.frame(familias)
generos <- as.data.frame(generos)

# Making sure every value is numeric
generos   <- apply(generos,c(1,2),function(x) 
  as.numeric(as.character(x)))
familias   <- apply(familias,c(1,2),function(x) 
  as.numeric(as.character(x)))

# Making a copy to work on, so "familias" and "generos" both remain untouched. 
# Also, exporting as .csv
countFamilias <- as.data.frame(familias)
countGeneros <- as.data.frame(generos)
write.csv(countFamilias, "D:/Laboratorio/Elisa/Microbiota2024/familias.csv",
          row.names=T)
write.csv(countGeneros, "D:/Laboratorio/Elisa/Microbiota2024/generos.csv",
          row.names=T)

# We need to remove features with a prevalence fo less than 10% (each feature 
# has to be in at least 10% of all samples). First we calculate, for each 
# feature, in how many samples they are not present.
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


anova_test(alpha_genus, Chao1 ~ Legend)                # p=0.01
tukey_hsd(alpha_genus, Chao1 ~ Legend)                 # Run-RunVel p=0.01
anova_test(alpha_genus, `Shannon Entropy` ~ Legend)    # p=0.041
tukey_hsd(alpha_genus, `Shannon Entropy` ~ Legend)     # Run-RunVel p=0.02
anova_test(alpha_genus, `Simpson Index` ~ Legend)      # p=0.056
tukey_hsd(alpha_genus, `Simpson Index` ~ Legend)       # n.s.

anova_test(alpha_families, Chao1 ~ Legend)             # n.s.
anova_test(alpha_families, `Shannon Entropy` ~ Legend) # n.s.
anova_test(alpha_families, `Simpson Index` ~ Legend)   # n.s.




plots <- list()
#estabas poniéndole las comparativas a pares en la alfa











plot_data <- alpha_genus[c(1:4)] %>%
  as.data.frame() %>%
  pivot_longer(!c("Legend"))



for (name in colnames(alpha_genus)[c(1:3)]) {
  # Filter data for the current name
  plot_data <- alpha_genus[c(1:4)] %>%
    as.data.frame() %>%
    pivot_longer(!c("Legend")) %>%
    filter(name == !!name)
  plot_data <- as.data.frame(plot_data)
  print(colnames(plot_data))
  maximo <- max(plot_data$value)
  minimo <- min(plot_data$value)
  rango <- maximo - minimo
  pes <- pairwise_t_test(plot_data, value ~ Legend, p.adjust.method = "fdr")
  print(pes)
  plot <-  ggboxplot(data = plot_data, x = "Legend", y = "value", 
                     fill = "Legend", outlier.shape = NA, 
                     palette = paletaEli, add = "dotplot", 
                     add.params = list(alpha = 0.5, dotsize = 0.75)) +
    stat_pvalue_manual(pes, tip.length = 0,
                       y.position = c(maximo + (rango * 0.05),
                                      maximo + (rango * 0.20),
                                      maximo + (rango * 0.10),
                                      maximo + (rango * 0.35),
                                      maximo + (rango * 0.28),
                                      maximo + (rango * 0.15))) +
    stat_compare_means(method = "anova", label.y = minimo - (rango * 0.10)) +
    ggtitle(paste(name)) +
    ylab("") + xlab("") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 12), axis.text.x = element_blank())
  
  # Store the plot in the list
  plots[[name]] <- plot
}

# Display the plots
wrap_plots(plots, ncol = 3, guides = "collect")






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

permanova_pairwise(dist_gen, grp = metadata$Group, permutations = 10000, 
                   method ="euclidean", padj = "holm")


############################################### Correlations


corGen <- generosCLR %>%
  t() %>%
  as.data.frame()  %>%
  add_column(MeanDiscriminationIndex = metadata$MeanDiscriminationIndex,
             .before = 1) 
rhoGen <- cor(corGen)
pcorGen <- cor_pmat(corGen)
signifCIGen <- which(as.data.frame(pcorGen)$MeanDiscriminationIndex <= 0.05)
signifCIGen <- rownames(as.data.frame(pcorGen))[c(signifCIGen)]

corFam <- familiasCLR %>%
  t() %>%
  as.data.frame()  %>%
  add_column(MeanDiscriminationIndex = metadata$MeanDiscriminationIndex,
             .before = 1) 
rhoFam <- cor(corFam)
pcorFam <- cor_pmat(corFam)
signifCIFam <- which(as.data.frame(pcorFam)$MeanDiscriminationIndex <= 0.05)
signifCIFam <- rownames(as.data.frame(pcorFam))[c(signifCIFam)]











familiasCLR %>%
  t() %>%
  as.data.frame()  %>%
  add_column(Group = metadata$Group, .before = 1) %>%
  add_column(MeanDiscriminationIndex = metadata$MeanDiscriminationIndex, .after = "Group") %>%
  pivot_longer(!c("Group", "MeanDiscriminationIndex"))  %>% 
  filter(name %in% signifCIFam)  %>% 
  ggscatter("value", "MeanDiscriminationIndex", 
            color = "Group",fill="lightgrey", palette=paletaEli) + 
  stat_cor(method = "pearson")+
  #facet_wrap(~name, scales="free_x") + 
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  geom_hline(yintercept = 0) + geom_smooth(method="lm", color="black") +
  ggtitle("Families")


generosCLR %>%
  t() %>%
  as.data.frame()  %>%
  add_column(Group = metadata$Group, .before = 1) %>%
  add_column(MeanDiscriminationIndex = metadata$MeanDiscriminationIndex, 
             .after = "Group") %>%
  pivot_longer(!c("Group", "MeanDiscriminationIndex"))  %>% 
  filter(name %in% signifCIGen)  %>% 
  ggscatter("value", "MeanDiscriminationIndex", 
            color = "Group",fill="lightgrey", palette=paletaEli) + 
  stat_cor(method = "pearson")+
  facet_wrap(~name, scales="free_x", ncol = 4) + geom_hline(yintercept = 0.2, 
                                                  linetype = "dashed") +
  geom_hline(yintercept = 0) + geom_smooth(method="lm", color="black")






# Correlation matrices
corGen <- genBH %>%
  t() %>%
  as.data.frame()  %>%
  add_column(MeanDiscriminationIndex = metadata$MeanDiscriminationIndex, 
             .before = 1) 

rhoGen <- cor(corGen)
pcorGen <- cor_pmat(corGen)
pcorGen <- as.data.frame(pcorGen)
rownames(pcorGen) <- rownames(rhoGen)
pcorGen$rowname <- NULL
ggcorrplot(rhoGen, p.mat = pcorGen, type = "lower", insig = "blank", lab=T, 
           colors = c("#a23203", "white", "#0373a2"))

# aquí las correlaciones de las variables cuyos anovas sean significativos
t(genBH)[,c(10,13,14,15)] %>% 
  as.data.frame()  %>%
  add_column(Group = metadata$Group, .before = 1) %>%
  add_column(MeanDiscriminationIndex = metadata$MeanDiscriminationIndex, 
             .after = "Group") %>%
  pivot_longer(!c("Group", "MeanDiscriminationIndex"))  %>% 
  ggscatter("value", "MeanDiscriminationIndex", 
            color = "Group",fill="lightgrey", palette=paletaEli) + 
  stat_cor(method = "pearson")+
  facet_wrap(~name, scales="free_x", ncol = 4) + 
  geom_hline(yintercept = 0.2, color="red") +
  geom_hline(yintercept = 0) + geom_smooth(method="lm", color="black")










# aquí las correlaciones de todos los géneros cuya correlación es significativa
# con el CI (independientemente de lo que saliese en su anova)
corGen <- generosCLR %>% t() %>% as.data.frame() %>% 
  add_column(MeanDiscriminationIndex = metadata$MeanDiscriminationIndex, 
             .before = 1)
rhoGen <- cor(corGen)
pcorGen <- cor_pmat(corGen)
pcorGen <- as.data.frame(pcorGen)
rownames(pcorGen) <- rownames(rhoGen)
pcorGen$rowname <- NULL
signif <- which(pcorGen$MeanDiscriminationIndex <= 0.05)

t(generosCLR)[,c(signif-1)] %>% 
  as.data.frame()  %>%
  add_column(Group = metadata$Group, .before = 1) %>%
  add_column(MeanDiscriminationIndex = metadata$MeanDiscriminationIndex, 
             .after = "Group") %>%
  pivot_longer(!c("Group", "MeanDiscriminationIndex"))  %>% 
  ggscatter("value", "MeanDiscriminationIndex", 
            color = "Group",fill="lightgrey", palette=paletaEli) + 
  stat_cor(method = "pearson")+
  facet_wrap(~name, scales="free_x", ncol = 4) + 
  geom_hline(yintercept = 0.2, color="red") +
  geom_hline(yintercept = 0) + geom_smooth(method="lm", color="black")














############################################### Differential abundance


genus.glm = fw_glm(x             = generosCLR,
                   f             = ~ Group, 
                   metadata      = metadata, 
                   adjust.method = "BH", 
                   order         = "atc")
family.glm = fw_glm(x             = familiasCLR,
                    f             = ~ Group, 
                    metadata      = metadata, 
                    adjust.method = "BH", 
                    order         = "atc")
glimpse(genus.glm)

hist(genus.glm$`tukeys.Run-Sed p adj`, xlim = c(0, 1), breaks = 20)
hist(genus.glm$`coefs.GroupRun Pr(>|t|)`, xlim = c(0, 1), breaks = 20)

signifRunSed <- which(genus.glm$`tukeys.Run-Sed p adj` <= 0.05)
signifRTRun <- which(genus.glm$`tukeys.RunTime-Run p adj` <= 0.05)
signifRVRun <- which(genus.glm$`tukeys.RunVel-Run p adj` <= 0.05)
signifRVRT <- which(genus.glm$`tukeys.RunVel-RunTime p adj`<= 0.05)
signifRTSed <- which(genus.glm$`tukeys.RunTime-Sed p adj` <= 0.05)
signifRVSed <- which(genus.glm$`tukeys.RunVel-Sed p adj` <= 0.05)
signifANOVAs <- which(genus.glm$`anovas.Group Pr(>F).BH` <= 0.05)


signifRunSed <- genus.glm$feature[c(signifRunSed)]
signifRTRun <- genus.glm$feature[c(signifRTRun)]
signifRVRun <- genus.glm$feature[c(signifRVRun)]
signifRVRT <- genus.glm$feature[c(signifRVRT)]
signifRTSed <- genus.glm$feature[c(signifRTSed)]
signifRVSed <- genus.glm$feature[c(signifRVSed)]



genBH <- generosCLR[genus.glm[genus.glm$`anovas.Group Pr(>F).BH` < 0.1,
                              "feature"],]




genBH %>%
  t() %>%
  as.data.frame()  %>%
  add_column(Group = metadata$Group) %>%
  pivot_longer(!c("Group"))  %>%
  ggplot(aes(x     = Group, 
             y     = value, 
             fill  = Group, 
             group = Group)) + 
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_beeswarm(size = 2, cex = 3) + 
  facet_wrap(~name, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = paletaEli) + 
  stat_compare_means(method="anova") +
  ggtitle("Genus") +
  ylab("") + xlab("") + theme_bw() + theme(text = element_text(size = 12)) 


facet_names <- unique(genBH$name)


plots <- list()

# Loop through each name and create a plot
for (name in facet_names) {
  # Filter data for the current name
  plot_data <- genBH %>%
    t() %>%
    as.data.frame() %>%
    add_column(Group = metadata$Group) %>%
    pivot_longer(!c("Group")) %>%
    filter(name == !!name)
  maximo <- max(plot_data$value)
  minimo <- min(plot_data$value)
  rango <- maximo - minimo
  pes <- pairwise_t_test(plot_data, value ~ Group, p.adjust.method = "BH")
  plot <-  ggboxplot(data = plot_data, x = "Group", y = "value", fill = "Group", 
                     outlier.shape = NA, palette = paletaEli, add = "dotplot", 
                     add.params = list(alpha = 0.5, dotsize = 1)) +
    stat_pvalue_manual(pes, tip.length = 0,
                       y.position = c(maximo + (rango * 0.05),
                                      maximo + (rango * 0.20),
                                      maximo + (rango * 0.10),
                                      maximo + (rango * 0.35),
                                      maximo + (rango * 0.28),
                                      maximo + (rango * 0.15))) +
    stat_compare_means(method = "anova", label.y = minimo - (rango * 0.10)) +
    ggtitle(paste(name)) +
    ylab("") + xlab("") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank()) + 
    theme(text = element_text(size = 12), axis.text.x = element_blank())

  # Store the plot in the list
  plots[[name]] <- plot
}

# Display the plots
wrap_plots(plots, ncol = 3, guides = "collect")


























significant <- c()
todoConTodo <- function(df) {
  tdf <- as.data.frame(t(df))
  tdf$Group <- metadata$Group
  significant <- c()
  for(i in colnames(tdf)[1:(ncol(tdf) - 1)]) {
    formula <- reformulate("Group", response = paste0("`", i, "`"))
    model <- aov(formula, data = tdf)
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
  geom_beeswarm(size = 2, cex = 3) + 
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
  geom_beeswarm(size = 2, cex = 3) + 
  facet_wrap(~name, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = paletaEli) + 
  stat_compare_means(method="anova") +
  ggtitle("Genus") +
  ylab("") + xlab("") + theme_bw() + theme(text = element_text(size = 12))







psG <- tax_glom(ps, taxrank = "Genus")

psSedRun <- phyloseq::subset_samples(
  psG,
  Group %in% c("Sed", "Run")
)

psSedRV <- phyloseq::subset_samples(
  psG,
  Group %in% c("Sed", "RunVel")
)

psSedRT <- phyloseq::subset_samples(
  psG,
  Group %in% c("Sed", "RunTime")
)

psRunRT <- phyloseq::subset_samples(
  psG,
  Group %in% c("Run", "RunTime")
)

psRunRV <- phyloseq::subset_samples(
  psG,
  Group %in% c("Run", "RunVel")
)

psRTRV <- phyloseq::subset_samples(
  psG,
  Group %in% c("RunTime", "RunVel")
)

aldexSedRun <- run_aldex(psSedRun, group="Group")






upsetG <- get_upset(psG, factorNames="Group")
upset(upsetG,  order.by="freq", empty.intersections="on")
upsetF <- get_upset(psF, factorNames="Group")
upset(upsetF,  order.by="freq", empty.intersections="on")


# ya que colnames(otu_table(psG)) == rownames(taxonomiaG) es todo TRUE,
unicos <- upsetG[rowSums(upsetG) == 1,]
unicosRun <- unicos[which(unicos$Run == 1),]
unicosRunVel <- unicos[which(unicos$RunVel == 1),]
unicosRunTime <- unicos[which(unicos$RunTime == 1),]
unicosSed <- unicos[which(unicos$Sed == 1),]

unicosRun <- which(rownames(taxonomiaG) %in% rownames(unicosRun))
unicosRun<- taxonomiaG$Genus[c(unicosRun)]

unicosRunTime <- which(rownames(taxonomiaG) %in% rownames(unicosRunTime))
unicosRunTime<- taxonomiaG$Genus[c(unicosRunTime)]

unicosRunVel <- which(rownames(taxonomiaG) %in% rownames(unicosRunVel))
unicosRunVel <- taxonomiaG$Genus[c(unicosRunVel)]

unicosSed <- which(rownames(taxonomiaG) %in% rownames(unicosSed))
unicosSed <- taxonomiaG$Genus[c(unicosSed)]








