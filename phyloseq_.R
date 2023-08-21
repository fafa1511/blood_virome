# 22/06

library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library(phyloseq.extended)
library(metagenomeSeq)
library(tidyr)
library(vegan)
library(gplots)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(scales)

library(DT)
library(gridExtra)
library(cowplot)
library(forcats)



setwd("~/work/work/dico/res/new_data")
# read data
abundance_mat = read.table("presence_abundance.csv", header = TRUE, sep = "\t")
tax_mat = read.table("taxo_bis.csv", header = TRUE, sep = "\t")
sample_mat = read.csv("metadata.csv", header = TRUE, sep = ";")

# remove ControlNeg3 NegControl2 plasmaBL PM stoolIT MF samples
abundance_mat[c(2, 3, 5, 21, 44, 47)]
abundance_mat <- abundance_mat[ , -c(2, 3, 5, 21, 44, 47)]   


# define row names
abundance_mat <- abundance_mat %>%
  tibble::column_to_rownames("X") 

tax_mat <- tax_mat %>%
  tibble::column_to_rownames("X") 

sample_mat <- sample_mat %>%
  tibble::column_to_rownames("X")


# remove controlNeg contigs
abundance_mat[c(15, 23508, 28227, 31474),]
abundance_mat <- abundance_mat[-c(15, 23508, 28227, 31474),]

tax_mat[c(15, 23508, 28227, 31474),]
tax_mat <- tax_mat[-c(15, 23508, 28227, 31474),]

# remove ControlNeg3 NegControl2 plasmaBL PM stoolIT MF samples
sample_mat[c(1, 2, 5, 30, 48, 54),]
sample_mat <- sample_mat[-c(1, 2, 5, 30, 48, 54) ,]   


# Transform into matrixes
abundance_mat <- as.matrix(abundance_mat)
tax_mat <- as.matrix(tax_mat)

# remet les spikes present
#tax_mat[15,1]
#tax_mat[15,1] <- "non"

# tax_mat[23508,1] #[1] "SPAdes_SampleLDS_NODE147"
# tax_mat[23508,1] <- "non"
# 
# tax_mat[28227,1] #[1] "SPAdes_SampleMM_NODE668"
# tax_mat[28227,1] <- "non"
# 
# tax_mat[31474,1] #[1] "SPAdes_SampleN_NODE681"
# tax_mat[31474,1] <- "non"



# Transform to phyloseq objects
abundance = otu_table(abundance_mat, taxa_are_rows = TRUE)
contig_info = tax_table(tax_mat)
sample_info = sample_data(sample_mat)
ab <- phyloseq(abundance,contig_info,sample_info)
ab
min(sample_sums(ab))
colnames(abundance)
rownames(sample_info)
# check
all(colnames(abundance) %in% rownames(sample_info))

all(rownames(abundance) %in% rownames(contig_info))



sample_sums(ab)
# recalcul ab
ab_re = transform_sample_counts(ab, function(x) x / sum(x))
# check
sample_sums(ab_re)



plot_bar(ab_re, fill = "contig_contaminant") + 
  geom_bar(aes(color=contig_contaminant, fill=contig_contaminant), stat="identity", position="stack")




######### Keep only decontaminant contig
contig_decontaminant <- subset_taxa(ab_re, contig_contaminant %in% c("non"))
contig_decontaminant

otu_table(contig_contaminant)[,"Stool_O"]
#s_decontaminant = sample_sums(contig_decontaminant)
#sum(s_decontaminant)/54
# check
otu_table(contig_decontaminant)["SPAdes_SampleA_NODE3", "StoolST"] # 0.008462721
sample_sums(contig_decontaminant) 


########### recalcul abundance
ab_decontaminant = transform_sample_counts(contig_decontaminant, function(x) x / sum(x))
# check
sample_sums(ab_decontaminant)

otu_table(ab_decontaminant)["SPAdes_SampleA_NODE3", "StoolST"] # 0.01054056



plot_bar(ab_decontaminant, fill = "contig_viral") + 
  geom_bar(aes(color=contig_viral, fill=contig_viral), stat="identity", position="stack")




# Keep only viral contig
contig_decontaminant_viral <- subset_taxa(ab_decontaminant, contig_viral %in% c("oui"))

# check
sample_sums(contig_decontaminant_viral)
plot_bar(contig_decontaminant_viral, fill = "contig_viral")

# recalcul abundance
ab_viral = transform_sample_counts(contig_decontaminant_viral, function(x) x / sum(x))


# check
sample_sums(ab_viral)

plot_bar(ab_viral, fill = "contig_viral")


taxa_sums(ab_viral)

# Heatmap
help(facet_grid)
p3 <- plot_heatmap(ab_viral,low = "yellow", high = "red", sample.order = "Statut", trans = log_trans(4)) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12, angle = -30, hjust = 0, vjust=0.5)) +
  facet_grid(. ~ Statut, scales='free') +
  theme(strip.text.x = element_text(size=16))
p3$scales$scales[[1]]$name <- "Viral contigs"
p3$scales$scales[[2]]$name <- ""
p3


# select healthy or crohn sample
plasma <- subset_samples(ab_viral, sampleStatut=="Plasma")

stool <- subset_samples(ab_viral, sampleStatut=="Stool")

plasma
stool

###############################
# Beta-diversity
dist_bc <- phyloseq::distance(ab_viral, method = "bray")
dist_bc

dist_jcc <- phyloseq::distance(ab_viral, method = "jaccard")
dist_jcc

# on BC 
set.seed(1115)
help(adonis2)
adonis_bc <- adonis2(formula = dist_bc~Statut, data = data.frame(sample_data(ab_viral)))
adonis_bc


res.adonis.bc <- paste("Permanova, ", "R2", " = ", round(adonis_bc$R2[1],3), ", p =", adonis_bc$`Pr(>F)`[1],"*")

# on JCC
set.seed(1115)
adonis_jcc <- adonis2(formula = dist_jcc~Statut, data = data.frame(sample_data(ab_viral)))
adonis_jcc

res.adonis.jcc <- paste("Permanova, ", "R2", " = ", round(adonis_jcc$R2[1],3), ", p =", adonis_jcc$`Pr(>F)`[1],"*")

#PCoA with BC
pal_step <- c("lightpink2","steelblue2")
ord_PCoA_bc <- ordinate(ab_viral, method = "PCoA", distance = dist_bc)
p1 <- plot_ordination(ab_viral, ord_PCoA_bc, color = "Statut") +
  scale_color_manual(values=pal_step) +
  labs(subtitle = res.adonis.bc) +
  theme_bw() +
  ggtitle("BRAY-CURTIS")
p1


#PCoA with Jaccard
ord_PCoA_jcc <- ordinate(ab_viral, method = "PCoA", distance = dist_jcc)
p2 <- plot_ordination(ab_viral, ord_PCoA_jcc, color = "Statut") +
  scale_color_manual(values=pal_step) +
  labs(subtitle = res.adonis.jcc) +
  theme_bw() +
  ggtitle("JACCARD") 
p2


# Alpha-diversity
alpha.diversity<-estimate_richness(ab_viral, measures = "Shannon")
alpha.diversity

res.kruskal <- alpha.diversity %>% kruskal_test(Shannon ~ sample_data(ab_viral)$Statut)
res.kruskal

res.kruskal.short <- paste("Kruskal-Wallis test, Chi2 = ", round(res.kruskal$statistic[1],3), ", p = ", round(res.kruskal$p[1],3))
p4 <- plot_richness(ab_viral, color = "Statut", measures = "Shannon", x="Statut") +
  geom_boxplot() +
  labs(subtitle = res.kruskal.short) +
  scale_color_manual(values=pal_step) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  ggtitle("SHANNON INDEX") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) 
p4




############################ ANALISE ECOLOGIQUE ############################


##############################################################################
#################### PLASMA
# Beta-diversity

dist_bc_plasma <- phyloseq::distance(plasma, method = "bray")


dist_jcc_plasma <- phyloseq::distance(plasma, method = "jaccard")


# on BC 
set.seed(1115)

adonis_bc_plasma <- adonis2(formula = dist_bc_plasma~Statut, data = data.frame(sample_data(plasma)))
adonis_bc_plasma

res.adonis.bc.plasma <- paste("Permanova, ", "R2", " = ", round(adonis_bc_plasma$R2[1],3), ", p =", adonis_bc_plasma$`Pr(>F)`[1],"***")

# on JCC
set.seed(1115)
adonis_jcc_plasma <- adonis2(formula = dist_jcc_plasma~Statut, data = data.frame(sample_data(plasma)))
adonis_jcc_plasma

res.adonis.jcc.plasma <- paste("Permanova, ", "R2", " = ", round(adonis_jcc_plasma$R2[1],3), ", p =", adonis_jcc_plasma$`Pr(>F)`[1],"***")

#PCoA with BC
ord_PCoA_bc_plasma <- ordinate(plasma, method = "PCoA", distance = dist_bc_plasma)
p1_plasma <- plot_ordination(plasma, ord_PCoA_bc_plasma, color = "Statut") +
  scale_color_manual(values=pal_step) +
  labs(subtitle = res.adonis.bc.plasma) +
  theme_bw() +
  ggtitle("BRAY-CURTIS")
p1_plasma


#PCoA with Jaccard
ord_PCoA_jcc_plasma <- ordinate(plasma, method = "PCoA", distance = dist_jcc_plasma)
p2_plasma <- plot_ordination(plasma, ord_PCoA_jcc_plasma, color = "Statut") +
  scale_color_manual(values=pal_step) +
  labs(subtitle = res.adonis.jcc.plasma) +
  theme_bw() +
  ggtitle("JACCARD") 
p2_plasma


# Alpha-diversity
alpha.diversity.plasma <-estimate_richness(plasma, measures = "Shannon")
alpha.diversity.plasma

res.kruskal.plasma <- alpha.diversity.plasma %>% kruskal_test(Shannon ~ sample_data(plasma)$Statut)
res.kruskal.plasma

res.kruskal.plasma <- paste("Kruskal-Wallis test, Chi2 = ", round(res.kruskal.plasma$statistic[1],3), ", p = ", round(res.kruskal.plasma$p[1],3))
p4_plasma <- plot_richness(plasma, color = "Statut", measures = "Shannon", x="Statut") +
  geom_boxplot() +
  labs(subtitle = res.kruskal.plasma) +
  scale_color_manual(values=pal_step) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  ggtitle("SHANNON INDEX") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) 
p4_plasma



#################################################################
########## STOOL
dist_bc_stool <- phyloseq::distance(stool, method = "bray")


dist_jcc_stool <- phyloseq::distance(stool, method = "jaccard")


# on BC 
set.seed(1115)

adonis_bc_stool <- adonis2(formula = dist_bc_stool~Statut, data = data.frame(sample_data(stool)))
adonis_bc_stool

res.adonis.bc.stool <- paste("Permanova, ", "R2", " = ", round(adonis_bc_stool$R2[1],3), ", p =", adonis_bc_stool$`Pr(>F)`[1])

# on JCC
set.seed(1115)
adonis_jcc_stool <- adonis2(formula = dist_jcc_stool~Statut, data = data.frame(sample_data(stool)))
adonis_jcc_stool

res.adonis.jcc.stool <- paste("Permanova, ", "R2", " = ", round(adonis_jcc_stool$R2[1],3), ", p =", adonis_jcc_stool$`Pr(>F)`[1])

#PCoA with BC
ord_PCoA_bc_stool <- ordinate(stool, method = "PCoA", distance = dist_bc_stool)
p1_stool <- plot_ordination(stool, ord_PCoA_bc_stool, color = "Statut") +
  scale_color_manual(values=pal_step) +
  labs(subtitle = res.adonis.bc.stool) +
  theme_bw() +
  ggtitle("BRAY-CURTIS")
p1_stool


#PCoA with Jaccard
ord_PCoA_jcc_stool <- ordinate(stool, method = "PCoA", distance = dist_jcc_stool)
p2_stool <- plot_ordination(stool, ord_PCoA_jcc_stool, color = "Statut") +
  scale_color_manual(values=pal_step) +
  labs(subtitle = res.adonis.jcc.stool) +
  theme_bw() +
  ggtitle("JACCARD") 
p2_stool


# Alpha-diversity
alpha.diversity.stool <-estimate_richness(stool, measures = "Shannon")
alpha.diversity.stool

res.kruskal.stool <- alpha.diversity.stool %>% kruskal_test(Shannon ~ sample_data(stool)$Statut)
res.kruskal.stool

res.kruskal.stool <- paste("Kruskal-Wallis test, Chi2 = ", round(res.kruskal.stool$statistic[1],3), ", p = ", round(res.kruskal.stool$p[1],3))
p4_stool <- plot_richness(stool, color = "Statut", measures = "Shannon", x="Statut") +
  geom_boxplot() +
  labs(subtitle = res.kruskal.stool) +
  scale_color_manual(values=pal_step) + 
  theme_bw() + 
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  ggtitle("SHANNON INDEX") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) 
p4_stool

############################ spike

# read data for spike analyze
setwd("~/work/work/dico/res/new_data")

abundance_mat = read.table("presence_abundance.csv", header = TRUE, sep = "\t")
tax_mat = read.table("taxo_bis.csv", header = TRUE, sep = "\t")
sample_mat = read.csv("metadata.csv", header = TRUE, sep = ";")

# remove ControlNeg3 NegControl2 plasmaBL PM stoolIT MF samples
abundance_mat[c(2, 3, 5, 21, 44, 47)]
abundance_mat <- abundance_mat[ , -c(2, 3, 5, 21, 44, 47)]   


# define row names
abundance_mat <- abundance_mat %>%
  tibble::column_to_rownames("X") 

tax_mat <- tax_mat %>%
  tibble::column_to_rownames("X") 

sample_mat <- sample_mat %>%
  tibble::column_to_rownames("X")


# remove ControlNeg3 NegControl2 plasmaBL PM stoolIT MF samples
sample_mat[c(1, 2, 5, 30, 48, 54),]
sample_mat <- sample_mat[-c(1, 2, 5, 30, 48, 54) ,]   


# Transform into matrixes
abundance_mat <- as.matrix(abundance_mat)
tax_mat <- as.matrix(tax_mat)

# remet les spikes present
tax_mat[15,1]
tax_mat[15,1] <- "non"

tax_mat[23508,1] #[1] "SPAdes_SampleLDS_NODE147"
tax_mat[23508,1] <- "non"

tax_mat[28227,1] #[1] "SPAdes_SampleMM_NODE668"
tax_mat[28227,1] <- "non"

tax_mat[31474,1] #[1] "SPAdes_SampleN_NODE681"
tax_mat[31474,1] <- "non"



# Transform to phyloseq objects
abundance = otu_table(abundance_mat, taxa_are_rows = TRUE)
contig_info = tax_table(tax_mat)
sample_info = sample_data(sample_mat)
ab <- phyloseq(abundance,contig_info,sample_info)
ab

# check
all(colnames(abundance) %in% rownames(sample_info))

all(rownames(abundance) %in% rownames(contig_info))



########## Keep only decontaminant contig ##########
contig_decontaminant <- subset_taxa(ab, contig_contaminant %in% c("non"))
# check
sample_sums(contig_decontaminant)


# recalcul abundance
ab_decontaminant = transform_sample_counts(contig_decontaminant, function(x) x / sum(x))
# check
sample_sums(ab_decontaminant)


############# Keep only viral contig #############
contig_decontaminant_viral <- subset_taxa(ab_decontaminant, contig_viral %in% c("oui"))
# check
sample_sums(contig_decontaminant_viral)

# recalcul abundance
ab_viral = transform_sample_counts(contig_decontaminant_viral, function(x) x / sum(x))
# check
sample_sums(ab_viral)


### SPP1
ab_viral
sample_names(ab_viral)

SPP1 = c()
for (i in 1:length(sample_names(ab_viral))){
  SPP1[i] = otu_table(ab_viral)["SPAdes_SampleA_NODE47",sample_names(ab_viral)[i]]
}

nom = c()
for (i in 1:length(sample_names(ab_viral))){
  nom[i] = sample_names(ab_viral)[i]
  print(nom[i])
}

nom
data = rbind(nom, SPP1)
data

### M13
# SPAdes_SampleLDS_NODE147
tax_table(ab_viral)["SPAdes_SampleLDS_NODE147"]
otu_table(ab)["SPAdes_SampleLDS_NODE147"]

M13_147 = c()
for (i in 1:length(sample_names(ab_viral))){
  M13_147[i] = otu_table(ab_viral)["SPAdes_SampleLDS_NODE147",sample_names(ab_viral)[i]]
}

nom_147 = c()
for (i in 1:length(sample_names(ab_viral))){
  nom_147[i] = sample_names(ab_viral)[i]
}


data_147 = rbind(nom_147, M13_147)
data_147

# SPAdes_SampleMM_NODE668
tax_table(ab_viral)["SPAdes_SampleMM_NODE668"]
otu_table(ab)["SPAdes_SampleMM_NODE668"]


M13_668 = c()
for (i in 1:length(sample_names(ab_viral))){
  M13_668[i] = otu_table(ab_viral)["SPAdes_SampleMM_NODE668",sample_names(ab_viral)[i]]
}

nom_668 = c()
for (i in 1:length(sample_names(ab_viral))){
  nom_668[i] = sample_names(ab_viral)[i]
}


data_668 = rbind(nom_668, M13_668)
data_668

#SPAdes_SampleN_NODE681

tax_table(ab_viral)["SPAdes_SampleN_NODE681"]
otu_table(ab)["SPAdes_SampleN_NODE681"]


M13_681 = c()
for (i in 1:length(sample_names(ab_viral))){
  M13_681[i] = otu_table(ab_viral)["SPAdes_SampleN_NODE681",sample_names(ab_viral)[i]]
}

nom_681 = c()
for (i in 1:length(sample_names(ab_viral))){
  nom_681[i] = sample_names(ab_viral)[i]
}


data_681 = rbind(nom_681, M13_681)
data_681

#################### 
otu_table(ab_viral)["SPAdes_SampleMB_NODE234", "StoolBR"]
