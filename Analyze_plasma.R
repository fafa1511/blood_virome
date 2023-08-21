# 16/05

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
library(pheatmap)

library(DT)
library(gridExtra)
library(cowplot)
library(scales)
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


# Transform to phyloseq objects
abundance = otu_table(abundance_mat, taxa_are_rows = TRUE)
contig_info = tax_table(tax_mat)
sample_info = sample_data(sample_mat)
ab <- phyloseq(abundance,contig_info,sample_info)
ab

# check
all(colnames(abundance) %in% rownames(sample_info))

all(rownames(abundance) %in% rownames(contig_info))



sample_sums(ab)
# recalcul ab
ab_re = transform_sample_counts(ab, function(x) x / sum(x))
# check
sample_sums(ab_re)

########################## decontaminant ##########################
contig_decontaminant <- subset_taxa(ab_re, contig_contaminant %in% c("non"))
sample_sums(contig_decontaminant)
# recalcul abundance
ab_decontaminant = transform_sample_counts(contig_decontaminant, function(x) x / sum(x))

# check
sample_sums(ab_decontaminant)


###################### Keep only viral contig ######################
contig_decontaminant_viral <- subset_taxa(ab_decontaminant, contig_viral %in% c("oui"))
sample_sums(contig_decontaminant_viral)
# recalcul abundance
ab_viral = transform_sample_counts(contig_decontaminant_viral, function(x) x / sum(x))
# check
sample_sums(ab_viral)
ab_viral



# select only plasma samples

plasma_sample <- subset_samples(ab_viral, sampleStatut=="Plasma")
plasma_sample

# abundance for each contig
plot_bar(plasma_sample, fill = "contig_viral")


# which contig is present
# MM  
for (i in 1:length(rownames(otu_table(plasma_sample)))){
  if (otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaMM"] != 0) {
    print(rownames(otu_table(plasma_sample))[i])
  }
  
}


# BK 
for (i in 1:length(rownames(otu_table(plasma_sample)))){
  if (otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaBK"] != 0) {
    print(otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaBK"])
  }
  
}


otu_table(plasma_sample)

for (i in 1:length(rownames(otu_table(plasma_sample)))){
  if (otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"Plasma_B_R1R2"] > 0.1) {
    print(otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"Plasma_B_R1R2"])
  }
  
}

for (i in 1:length(rownames(otu_table(plasma_sample)))){
  if (otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"Plasma_E_R1R2"] > 0.1) {
    print(otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"Plasma_E_R1R2"])
  }
  
}

for (i in 1:length(rownames(otu_table(plasma_sample)))){
  if (otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"Plasma_K_R1R2"] > 0.1) {
    print(otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"Plasma_K_R1R2"])
  }
  
}


for (i in 1:length(rownames(otu_table(plasma_sample)))){
  if (otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"Plasma_L_R1R2"] > 0.1) {
    print(otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"Plasma_L_R1R2"])
  }
  
}

for (i in 1:length(rownames(otu_table(plasma_sample)))){
  if (otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaCMD"] > 0.1) {
    print(otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaCMD"])
  }
  
}


for (i in 1:length(rownames(otu_table(plasma_sample)))){
  if (otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaLDS"] > 0.1) {
    print(otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaLDS"])
  }
  
}



for (i in 1:length(rownames(otu_table(plasma_sample)))){
  if (otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaMM"] > 0.1) {
    print(otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaMM"])
  }
  
}

for (i in 1:length(rownames(otu_table(plasma_sample)))){
  if (otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaNR"] > 0.1) {
    print(otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaNR"])
  }
  
}

for (i in 1:length(rownames(otu_table(plasma_sample)))){
  if (otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaBK"] > 0.1) {
    print(otu_table(plasma_sample)[rownames(otu_table(plasma_sample))[i],"PlasmaBK"])
  }
  
}



# Heatmap
setwd("~/work/work/dico/res/new_data")
# read data
viral_presence_mat = read.table("viral_presence_plasma.csv", header = TRUE, sep = "\t")
tax_mat = read.table("taxo_bis.csv", header = TRUE, sep = "\t")
sample_viral_presece_mat = read.csv("metadata_viral_presence.csv", header = TRUE, sep = ";")
# remove PlasmaBL PM samples
viral_presence_mat[c(3, 19)]
viral_presence_mat <- viral_presence_mat[ , -c(3, 19)]

# define row names
viral_presence_mat <- viral_presence_mat %>%
  tibble::column_to_rownames("X")

tax_mat <- tax_mat %>%
  tibble::column_to_rownames("X")

sample_viral_presece_mat <- sample_viral_presece_mat %>%
  tibble::column_to_rownames("X")


# remove controlNeg contigs
viral_presence_mat[c(15, 23508, 28227, 31474),]
viral_presence_mat <- viral_presence_mat[-c(15, 23508, 28227, 31474),]

tax_mat[c(15, 23508, 28227, 31474),]
tax_mat <- tax_mat[-c(15, 23508, 28227, 31474),]

# remove PlasmaBL PM
sample_viral_presece_mat[c(3, 28),]
sample_viral_presece_mat <- sample_viral_presece_mat[-c(3, 28) ,]


# Transform into matrixes
viral_presence_mat <- as.matrix(viral_presence_mat)
tax_mat <- as.matrix(tax_mat)




# Transform to phyloseq objects
presence = otu_table(viral_presence_mat, taxa_are_rows = TRUE)
contig_info = tax_table(tax_mat)
sample_info = sample_data(sample_viral_presece_mat)
presence <- phyloseq(presence,contig_info,sample_info)

# check
all(colnames(presence) %in% rownames(sample_info))

all(rownames(presence) %in% rownames(contig_info))

presence

# remove row sum = 0
ind = c()
for (i in 1:length(rownames(otu_table(presence)))){
  if (sum(otu_table(presence)[i]) == 0) {
    ind[i] = i
  }
}


# count NA
sum(is.na(ind)) 
# [1] 1794


# remove NA : omit NA values from vector
ind_bis <- as.numeric(na.omit(ind))
ind_bis
length(ind_bis) 
# [1] 84187

# remove 0 lines from phyloseq
presence 
# [ 85981 taxa and 27 samples ]
otu_table(presence)

otu_table(presence) <- otu_table(presence)[-ind_bis]

presence
# [ 1794 taxa and 27 samples ]



presence

# convert physeq obj to dataframe
df_presence = as.data.frame(as(otu_table(presence), "matrix"))

# have list of sum of each line
som = as.numeric(rowSums(df_presence))
som

# add som column to dataframe
df_presence$somme = som


# sort dataframe by column 'som'
df_presence
df_trier <- df_presence[order(df_presence[, 28], decreasing = TRUE),]
df_trier
df_som <- df_trier
df_som[28]
df_trier <- df_trier[,-28]

# count the number of contigs for each value of the sum
tab = rep(0, length(df_som[28]$somme))
table(df_som[28]$somme, tab)
list_col = as.numeric(table(df_som[28]$somme, tab))
list_col



annotation_row_bis = data.frame(
  Sommme_contig = factor(df_som[,28])
)
rownames(annotation_row_bis) = rownames(df_som[28])
annotation_row_bis

ann_colors_bis = list(
  Sommme_contig = c("1" = "#7570B3", "2" = "#8ECEE9", "3" = "#ADD39D", "4" = "#FFB888", "5" = "#FFB6F1", 
            "6" = "#FCE37A", "7" = "#FCE37A", "8" = "#FCE37A", "9" = "#E7298A", "10" = "#E7298A",
            "11" = "#E7298A", "12" = "#E7298A", "13" = "#E7298A", "14" = "#E7298A", "15" = "#E7298A",
            "16" = "#E7298A", "17" = "#E7298A")
)
ann_colors_bis


df_trier
df_som[,28]
rownames(df_som[28])


pheatmap(df_trier,
         color = c("black","orange"),
         cluster_cols = TRUE,# pour faire de réarrangement de colonnes
         clustering_distance_cols = "euclidean",
         cluster_rows = FALSE,# pour ne pas faire de réarrangement de lignes
         legend = FALSE, # Ne pas mettre la legende
         show_rownames = FALSE,# Ne pas mettre row names
         annotation_row = annotation_row_bis,
         annotation_legend = FALSE,
         annotation_colors = ann_colors_bis) 


# 7 contigs with high level of abundance
df = df_som[28]

# select all 995 lines (somme > 1)
rownames(df)[1:995]

"SPAdes_SampleB_NODE324"%in%rownames(df)[1:995]
"SPAdes_SampleK_NODE3"%in%rownames(df)[1:995]
"SPAdes_SampleCMD_NODE1737"%in%rownames(df)[1:995]
"SPAdes_SampleBK_NODE246"%in%rownames(df)[1:995]

# contig shared by more than 2 people
"SPAdes_SampleA_NODE1490"%in%rownames(df)[1:995]
"SPAdes_SampleL_NODE115"%in%rownames(df)[1:995]
"Cross_Assembly_MB_NODE21128"%in%rownames(df)[1:995]

rownames(df)[1:102]

sum(otu_table(presence)["SPAdes_SampleA_NODE1490"]) # B, 0.1784323 
# [1] 4

sum(otu_table(presence)["SPAdes_SampleL_NODE115"]) # L, 0.3003997
# [1] 10

sum(otu_table(presence)["Cross_Assembly_MB_NODE21128"]) # NR, 0.1230005
# [1] 14

# list of transfer contigs in shared contigs list (more than 10 people)
transfert = c("SPAdes_SampleA_NODE1177", "SPAdes_SamplePM_NODE1617", "SPAdes_SampleST_NODE1913", "SPAdes_SampleL_NODE30", "SPAdes_SampleL_NODE5485", "Cross_Assembly_MB_NODE3237", "Cross_Assembly_MB_NODE8120", "SPAdes_SampleL_NODE30", "Cross_Assembly_MB_NODE1712", "SPAdes_SampleBO_NODE6328", "SPAdes_SampleBM_NODE922", "SPAdes_SampleST_NODE702", "SPAdes_SampleBP_NODE2059", "SPAdes_SampleBP_NODE727", "Cross_Assembly_MB_NODE15411", "SPAdes_SampleB_NODE566")
for (i in 1:length(transfert)){
  print(transfert[i]%in%rownames(df)[1:102])
}
# [1] FALSE
# [1] FALSE
# [1] FALSE
# [1] FALSE
# [1] FALSE
# [1] FALSE
# [1] FALSE
# [1] FALSE
# [1] FALSE
# [1] TRUE  SPAdes_SampleBO_NODE6328
# [1] FALSE
# [1] FALSE
# [1] FALSE
# [1] FALSE
# [1] FALSE
# [1] FALSE

sum(otu_table(presence)["SPAdes_SampleBO_NODE6328"])  # [1] 11



