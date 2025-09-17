###############################################################################################
# Investigating the introduction route of the nonnative yellowbar angelfish Pomacanthus maculosus 
#in the southwestern Atlantic through otolith chemistry
###############################################################################################

#Authors: Johnatas Adelir-Alves, Patrick Derviche*, Claudio Oliveira, Daphne Spier, Humberto Gerum,
#Felippe Alexandre Daros, Beatriz Rochitti Boza, Marcelo Soeth

###################
#R script made by Patrick Derviche
###################


#: Summary
#: 1. Read dataset and clean
#: 2. Transform and summarize data (ratios, age groups)
#: 3. PERMANOVA for age-group effects
#: 4. PERMDISP excluding Age 0
#: 5. Pairwise Adonis post hoc tests
#: 6. PCA with visualization


###############################

# Working directory
setwd("C:/Users/patri/OneDrive/Documentos/Pommac")

# Load libraries
if(!require(pacman)) install.packages("pacman")
pacman::p_load("devtools","ggplot2","ggpubr","vegan","dplyr","corrplot","ggbiplot","scatterplot3d",
               "ggpmisc","FactoMineR","factoextra","GGally","MVN","cluster","ggdendro","gplots","NbClust")

###############################
# 1. Read dataset
###############################

# Clean R environment 
rm(list = ls())

# Import dataset
data <-read.csv2("data_POMMAC.csv")

# Select columns
data <- data[, c("Age","Mg24", "Sr87", "Ba138")]

# Ensure numeric types
data$Ba138 <- as.numeric(data$Ba138)
data$Sr87 <- as.numeric(data$Sr87)
data$Mg24 <- as.numeric(data$Mg24)
data$Age <- as.character(data$Age)

str(data)



###############################
# 2. Summary statistics
############################### 

# Scale values
data$Mg24 <- (data$Mg24) * 100000
data$Sr87 <- (data$Sr87) * 1000
data$Ba138 <- (data$Ba138) * 100000

# Summary table
data.frame(
  Element = c("Mg24", "Sr87", "Ba138"),
  Mean = c(round(mean(data$Mg24, na.rm = TRUE), 3),
           round(mean(data$Sr87, na.rm = TRUE), 3),
           round(mean(data$Ba138, na.rm = TRUE), 3)),
  SD = c(round(sd(data$Mg24, na.rm = TRUE), 3),
         round(sd(data$Sr87, na.rm = TRUE), 3),
         round(sd(data$Ba138, na.rm = TRUE), 3)),
  Min = c(round(min(data$Mg24, na.rm = TRUE), 3),
          round(min(data$Sr87, na.rm = TRUE), 3),
          round(min(data$Ba138, na.rm = TRUE), 3)),
  Max = c(round(max(data$Mg24, na.rm = TRUE), 3),
          round(max(data$Sr87, na.rm = TRUE), 3),
          round(max(data$Ba138, na.rm = TRUE), 3)))

#The Mg/Ca from 3.720 to 18.449 μg g−1 (5.700 ± 2.327 μg g−1; mean ± standard 
#deviation), the Sr/Ca ratios ranged from 3.976 to 10.882 mg g−1
#(6.933 ± 1.669 mg g−1), and the Ba/Ca from 0.469 to 7.140 μg g−1 (2.287 ± 1.691 mg g−1).

# Rename columns to ratios
data <- data %>%
  rename(
    `Mg:Ca` = Mg24,
    `Sr:Ca` = Sr87,
    `Ba:Ca` = Ba138 )

str(data)

# Add age groups
data_group <- data %>%
  mutate(`Age_group` = case_when(
    Age %in% 0:2 ~ "0-2",
    Age %in% 3:4 ~ "3-4",
    Age %in% 5:6 ~ "5-6",
    Age %in% 7:8 ~ "7-8",
    Age %in% 9:10 ~ "9-10",
    TRUE ~ NA_character_))

str(data_group)



###############################
# 3. PERMANOVA
###############################

# Select variables
summary(data_group)
var <- data_group[, c("Mg:Ca", "Sr:Ca", "Ba:Ca")]
summary(var)

# Apply log transformation
var.log <- log(var)

# Replace the original data with the transformed versions
data_group[, c("Mg:Ca", "Sr:Ca", "Ba:Ca")] <- var.log

# Perform adonis2
data.matrix <- as.matrix(var.log) 
data.Zmatrix <- scale(data.matrix)
matrix_dist <- vegdist(data.Zmatrix, method='euclidean')
permanova <- adonis2(matrix_dist ~ Age_group, data=data_group, permutations = 999, method="euclidean")

# Verify significant difference
permanova 

#The PERMANOVA indicated significant differences in otolith chemistry according to age (F = 23.976, p value < 0.001).



###############################
# 4. PERMDISP
############################### 

# Group data
data_group_without00 <- data %>%
  mutate(`Age_group` = case_when(
    Age %in% 1:2 ~ "1-2",
    Age %in% 3:4 ~ "3-4",
    Age %in% 5:6 ~ "5-6",
    Age %in% 7:8 ~ "7-8",
    Age %in% 9:10 ~ "9-10",
    TRUE ~ NA_character_))

data_group_without00 <- na.omit(data_group_without00)

summary(data_group_without00)
var_without00 <- data_group_without00[, c("Mg:Ca", "Sr:Ca", "Ba:Ca")]
summary(var_without00)

# Apply log transformation (natural log)
var.log_without00 <- log(var_without00)

# Replace the original data with the transformed versions
data_group_without00[, c("Mg:Ca", "Sr:Ca", "Ba:Ca")] <- var.log_without00


data.matrix_without00 <- as.matrix(var.log_without00) 
data.Zmatrix_without00 <- scale(data.matrix_without00)
matrix_dist_without00 <- vegdist(data.Zmatrix_without00, method='euclidean')
permanova_without00 <- adonis2(matrix_dist_without00 ~ Age_group, data=data_group_without00, permutations = 999, method="euclidean")

# Results
dispersion_without00 <- betadisper(matrix_dist_without00, group = data_group_without00$Age_group)
disp_test_without00 <- permutest(dispersion_without00, permutations = 999)
disp_test_without00

# Visualizations
plot(dispersion_without00)
boxplot(dispersion_without00, main = "Multivariate Dispersion (by Island)", ylab = "Distance to Centroid")


###############################
# 5. post hoc Pairwise Adonis
###############################

#Teste post hoc
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pairwise.adonis(matrix_dist, data_group$Age_group, sim.method = "euclidean",
                p.adjust.m = "bonferroni")

#The post hoc Adonis’ pairwise test indicated that the age class 0–2 years differed 
#significantly from 3–4, 5–6, and 7–8 years (p adjusted < 0.01), and showed marginal 
#differences compared to 9–10 years (p adjusted = 0.02). Conversely, the oldest group (9–10 years) 
#differed significantly only from the intermediate classes 3–4, 5–6, and 7–8 years (p adjusted < 0.01). 
#No significant differences were detected among the intermediate classes (3–4, 5–6, and 7–8 years; p adjusted ≥ 0.11).



###############################
# 6. PCA
###############################

# PCA
data$Age <- as.numeric(data$Age)

res.pca <- PCA(data[, 2:4], graph = FALSE,    
               scale.unit = TRUE,           
               ncp = 3)                    
res.pca

# Eigenvalues
round(get_eigenvalue(res.pca),3)  

# Groups
groups <- as.factor(data$Age[1:142])
dim(data)

# Figure
fviz_pca_biplot(res.pca, 
                col.ind = groups, 
                addEllipses = TRUE, 
                ellipse.type = "confidence",  ellipse.level = 0.95,
                legend.title = "Age",
                geom.ind = "point",
                pointsize = 1.5,
                repel = TRUE,
                alpha.var = 0.5,
                alpha.ind = 0.3,
                title = " ")
+  theme_light() 


###############################
# END
###############################
# Patrick Derviche