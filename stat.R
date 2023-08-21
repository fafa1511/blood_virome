# Analyze statistic of data

library("ggsignif")
library('ggpubr')
library("dplyr")

############################ data plasma ############################
# without PM/BL
nb_contig_p <- c(368,120,552,100,132,6,110,53,66,240,1,63,585,59,94,140,173,250,11,257,183,300,223,370,101,115,571)
sample_p <- c("A","B","C","D","E","CMD","IT","LDS","MB","MF","MM","NR","ST","K","L","M","N","O","BK","BM","BN","BO","BP","BQ","BR","BS", "BT")
statut_p <- c("Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy","Healthy","Healthy","Healthy","Healthy","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn")

mydata_p <- data.frame(sample = sample_p, statut = statut_p, nb_contig = nb_contig_p)
mydata_p

summary(mydata_p)
head(mydata_p)

# Calculate summary statistics
group_by(mydata_p, statut) %>%
  summarise(
    count = n(),
    mean = mean(nb_contig, na.rm = TRUE),
    sd = sd(nb_contig, na.rm = TRUE)
  )

# before doing the t-test, check the distribution is normal or not
ggqqplot(mydata_p$nb_contig)


# Shapiro tes
shapiro.test(nb_contig_p) # p-value = 0.002053
# Shapiro-Wilk normality test for Healthy's nb_contig
with(mydata_p, shapiro.test(nb_contig[statut == "Healthy"]))# p-value = 0.006867
# Shapiro-Wilk normality test for Crohn's nb_contig
with(mydata_p, shapiro.test(nb_contig[statut == "Crohn"])) # p-value = 0.209



# Ficher test for variance
res.ftest <- var.test(nb_contig ~ statut, data = mydata_p) 
res.ftest # p-value = 0.2844

# t-test
res <- t.test(nb_contig ~ statut, data = mydata_p, var.equal = TRUE)
res # p-value = 0.7753

# boxplot
my_comparisons <- list( c("Healthy", "Crohn") )

ggboxplot(mydata_p, x = "statut", y = "nb_contig",
          xlab = "statut des individus",
          ylab = "nombre de contig dans le sang",
          bxp.errorbar = TRUE,
          color = "statut", 
          palette =c("steelblue2", "lightpink2"),
          add = "jitter",
          shape = "statut",
          order = c("Healthy", "Crohn"))+
  ## add global test p-value
  stat_compare_means(method = "t.test", label.x = "Crohn", label.y=650)+
  ## add pairwise comparison p values
  geom_signif(comparisons = my_comparisons,
              map_signif_level = TRUE,
              y_position = c(607.5, 609))+
  theme(legend.position="none") 


############################ data stool ############################
# without IT/MF
nb_contig_s <- c(2094,660,2064,2189,1558,828,724,788,3186,666,1065,1091,196,2486,1136,1515,323,166,95,750,230,1091,1180,163,344,76,546)
sample_s <- c("A","B","C","D","E","CMD","LDS","MB","MM","NR","PM", "ST","K","L","M","N","O","BK","BL","BM","BN","BO","BP","BQ","BR","BS", "BT")
statut_s <- c("Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy","Healthy", "Healthy","Healthy","Healthy","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn")

mydata_s <- data.frame(sample = sample_s, statut = statut_s, nb_contig = nb_contig_s)
mydata_s

summary(mydata_s)

# Calculate summary statistics
group_by(mydata_s, statut) %>%
  summarise(
    count = n(),
    mean = mean(nb_contig, na.rm = TRUE),
    sd = sd(nb_contig, na.rm = TRUE)
  )



# Distribution
ggqqplot(mydata_s$nb_contig) # suit loi normale ou pas

# Shapiro-Wilk normality test for Healthy's nb_contig
with(mydata_s, shapiro.test(nb_contig[statut == "Healthy"]))# p-value = 0.04111
# Shapiro-Wilk normality test for Crohn's nb_contig
with(mydata_s, shapiro.test(nb_contig[statut == "Crohn"])) # p-value = 0.008438



# Ficher test for variance
res.ftest <- var.test(nb_contig ~ statut, data = mydata_s) 
res.ftest  # p-value = 0.5399

# Compute t-test
res <- t.test(nb_contig ~ statut, data = mydata_s, var.equal = TRUE)
res # p-value = 0.01856


# boxplot
my_comparisons <- list( c("Healthy", "Crohn") )

ggboxplot(mydata_s, x = "statut", y = "nb_contig", 
          xlab = "statut des individus",
          ylab = "nombre de contig dans la fèces ",
          bxp.errorbar = TRUE,
          color = "statut", 
          palette =c("steelblue2", "lightpink2"),
          add = "jitter",
          shape = "statut",
          order = c("Healthy", "Crohn"))+
  ## add global test p-value
  stat_compare_means(method = "t.test", label.x = "Crohn", label.y=3500)+
  ## add pairwise comparison p values
  geom_signif(comparisons = my_comparisons,
              map_signif_level = TRUE,
              y_position = c(3207.5, 3209))+
  theme(legend.position="none")



############################ data contig en commun ############################
# without IT/MF/PM/BK
nb_contig_c <- c(2,0,0,0,0,0,0,0,0,0,1,1,3,0,2,1,0,2,0,0,1,2,1,1,0)
sample_c <- c("A","B","C","D","E","CMD","LDS","MB","MM","NR", "ST","K","L","M","N","O","BK","BM","BN","BO","BP","BQ","BR","BS", "BT")
statut_c <- c("Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy","Healthy", "Healthy","Healthy","Healthy","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn")

mydata_c <- data.frame(sample = sample_c, statut = statut_c, nb_contig = nb_contig_c)
mydata_c

summary(mydata_c)

# Calculate summary statistics
group_by(mydata_c, statut) %>%
  summarise(
    count = n(),
    mean = mean(nb_contig, na.rm = TRUE),
    sd = sd(nb_contig, na.rm = TRUE)
  )



# Distribution
ggqqplot(mydata_c$nb_contig) 

# Shapiro-Wilk normality test for Healthy's nb_contig
with(mydata_c, shapiro.test(nb_contig[statut == "Healthy"]))# p-value = 1.687e-06
# Shapiro-Wilk normality test for Crohn's nb_contig
with(mydata_c, shapiro.test(nb_contig[statut == "Crohn"])) # p-value = 0.0318



# Ficher test for variance
res.ftest <- var.test(nb_contig ~ statut, data = mydata_c) 
res.ftest # p-value = 0.2152 -> variance not equal for two groups -> Wilcoxon test


# Wilcoxon test
res_wilcox <- wilcox.test(nb_contig ~ statut, data = mydata_c,
                          exact = FALSE)
res_wilcox # p-value = 0.03301


# boxplot
mydata_c$statut <- factor(mydata_c$statut, levels = c("Healthy", "Crohn"), ordered = TRUE)

new_plot <- ggplot(mydata_c, aes(x = statut, y = nb_contig, color = statut, shape = statut)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar") +
  geom_point(position = position_jitter(w = 0.1, h = 0), size = 2) +
  scale_color_manual(values = c("steelblue2", "lightpink2")) +
  xlab("statut des individus") +
  ylab("nombre de contig partage en commun") +
  stat_compare_means(method = "wilcox.test", label.x = "Crohn", label.y=3.3) +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = TRUE,
              y_position = c(3.1, 39),
              color = "black") +
  theme_classic() +
  theme(legend.position="none")
new_plot


############################ %/plasma ############################
# without IT/MF/PM/BK
pourcentage_contig_pp <- c(0.543,0,0,0,0,0,0,0,0,0,0.171,1.695,3.191,0,1.156,0.400,0,0.778,0,0,0.448,0.541,0.990,0.870,0)
sample_pp <- c("A","B","C","D","E","CMD","LDS","MB","MM","NR", "ST","K","L","M","N","O","BK","BM","BN","BO","BP","BQ","BR","BS", "BT")
statut_pp <- c("Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy","Healthy", "Healthy","Healthy","Healthy","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn")

mydata_pp <- data.frame(sample = sample_pp, statut = statut_pp, pourcentage_contig = pourcentage_contig_pp)
mydata_pp
summary(mydata_pp)

# Calculate summary statistics
group_by(mydata_pp, statut) %>%
  summarise(
    count = n(),
    mean = mean(pourcentage_contig, na.rm = TRUE),
    sd = sd(pourcentage_contig, na.rm = TRUE)
  )



# Distribution
ggqqplot(mydata_pp$pourcentage_contig) 

# Shapiro-Wilk normality test for Healthy's nb_contig
with(mydata_pp, shapiro.test(pourcentage_contig[statut == "Healthy"]))# p-value = 6.373e-07
# Shapiro-Wilk normality test for Crohn's nb_contig
with(mydata_pp, shapiro.test(pourcentage_contig[statut == "Crohn"])) # p-value = 0.004214



# Ficher test for variance
res.ftest <- var.test(pourcentage_contig ~ statut, data = mydata_pp)  
res.ftest # p-value = 8.18e-06 -> variance not equal for two groups -> Wilcoxon test

# Wilcoxon test
res_wilcox <- wilcox.test(pourcentage_contig ~ statut, data = mydata_pp,
                          exact = FALSE)
res_wilcox # p-value = 0.01347


# boxplot
my_comparisons <- list( c("Healthy", "Crohn") )

ggboxplot(mydata_pp, x = "statut", y = "pourcentage_contig", 
          xlab = "statut des individus",
          ylab = "% de contig en commun/contig tot dans le sang",
          bxp.errorbar = TRUE,
          color = "statut", 
          palette =c("steelblue2", "lightpink2"),
          add = "jitter",
          shape = "statut",
          order = c("Healthy", "Crohn"))+
  ## add global test p-value
  stat_compare_means(method = "wilcox.test", label.x = "Crohn", label.y=2.4)+
  ## add pairwise comparison p values
  geom_signif(comparisons = my_comparisons,
              map_signif_level = TRUE,
              y_position = c(2.2, 29))+
  theme(legend.position="none")


############################ %/stool ############################
# without IT/MF/PM/BK
pourcentage_contig_ps <- c(0.096,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.092,0.510,0.121,0.000,0.132,0.310,0.000,0.267,0.000,0.000,0.085,1.227,0.291,1.316,0)
sample_ps <- c("A","B","C","D","E","CMD","LDS","MB","MM","NR", "ST","K","L","M","N","O","BK","BM","BN","BO","BP","BQ","BR","BS", "BT")
statut_ps <- c("Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Healthy","Healthy", "Healthy","Healthy","Healthy","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn","Crohn")

mydata_ps <- data.frame(sample = sample_ps, statut = statut_ps, pourcentage_contig = pourcentage_contig_ps)
mydata_ps

summary(mydata_ps)

# Calculate summary statistics
group_by(mydata_ps, statut) %>%
  summarise(
    count = n(),
    mean = mean(pourcentage_contig, na.rm = TRUE),
    sd = sd(pourcentage_contig, na.rm = TRUE)
  )


# Distribution
ggqqplot(mydata_ps$pourcentage_contig)

# Shapiro-Wilk normality test for Healthy's nb_contig
with(mydata_ps, shapiro.test(pourcentage_contig[statut == "Healthy"]))# p-value = 1.174e-06
# Shapiro-Wilk normality test for Crohn's nb_contig
with(mydata_ps, shapiro.test(pourcentage_contig[statut == "Crohn"])) # p-value = 0.0004536



# Ficher test for variance
res.ftest <- var.test(pourcentage_contig ~ statut, data = mydata_ps)  
res.ftest # p-value = 4.287e-09 -> variance not equal for two groups -> Wilcoxon test


# Wilcoxon test
res_wilcox <- wilcox.test(pourcentage_contig ~ statut, data = mydata_ps,
                          exact = FALSE)
res_wilcox # p-value = 0.01136

# boxplot
my_comparisons <- list( c("Healthy", "Crohn") )

ggboxplot(mydata_ps, x = "statut", y = "pourcentage_contig", 
          xlab = "statut des individus",
          ylab = "% de contig en commun/contig tot dans la fèces ",
          bxp.errorbar = TRUE,
          color = "statut", 
          palette =c("steelblue2", "lightpink2"),
          add = "jitter",
          shape = "statut",
          order = c("Healthy", "Crohn"))+
  ## add global test p-value
  stat_compare_means(method = "wilcox.test", label.x = "Crohn", label.y=1.5)+
  ## add pairwise comparison p values
  geom_signif(comparisons = my_comparisons,
              map_signif_level = TRUE,
              y_position = c(1.4, 19))+
  theme(legend.position="none")