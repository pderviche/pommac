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
#: 7. Generalized Additive Model and 95% Bayesian confidence intervals


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

###
# Age classes
###

# PCA
res.pca <- PCA(data_group[, 2:4], graph = FALSE,    
               scale.unit = TRUE,           
               ncp = 3)                    
res.pca

# Eigenvalues
round(get_eigenvalue(res.pca),3)  

# Groups
groups <- as.factor(data_group$Age_group[1:142])
dim(data)

# Figure
fviz_pca_biplot(res.pca, 
                col.ind = groups, 
                addEllipses = TRUE, 
                ellipse.type = "confidence",  ellipse.level = 0.95,
                legend.title = "Age classes",
                geom.ind = "point",
                pointsize = 1.5,
                repel = TRUE,
                alpha.var = 0.5,
                alpha.ind = 0.3,
                title = " ")
+  theme_light() 

###
# Individual age
###

# PCA
data$Age <- as.numeric(data_group$Age)

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
p <- fviz_pca_biplot(res.pca, 
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

p

ggsave("Figure 3.png", p, width = 5, height = 5, dpi = 300, bg = "white")

###############################
# 7. Generalized Additive Model and 95% Bayesian confidence intervals
###############################

# Packages
library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# -----------------------------
# 1) Import dataset
# -----------------------------
data <- read.csv2("data_POMMAC.csv")

# Select columns and set types
data <- data[, c("Measure","Age","Mg24","Sr87","Ba138")]

data <- data %>%
  mutate(
    Ba138   = as.numeric(Ba138),
    Sr87    = as.numeric(Sr87),
    Mg24    = as.numeric(Mg24),
    Measure = round(as.numeric(Measure), 3),
    Age     = as.character(Age)
  )

# Convert to µmol·mol⁻¹
data <- data %>%
  mutate(
    Ba138 = Ba138 * 1000,
    Mg24  = Mg24  * 1000,
    Sr87  = Sr87  * 1000
  )

# -----------------------------
# 2) GAM models (positive data)
# -----------------------------
m_ba <- gam(Ba138 ~ s(Measure, k = 20), data = data,
            method = "REML", family = Gamma(link = "log"))
m_mg <- gam(Mg24  ~ s(Measure, k = 20), data = data,
            method = "REML", family = Gamma(link = "log"))
m_sr <- gam(Sr87  ~ s(Measure, k = 20), data = data,
            method = "REML", family = Gamma(link = "log"))

# Prediction grid
newx <- data.frame(
  Measure = seq(min(data$Measure, na.rm = TRUE),
                max(data$Measure, na.rm = TRUE),
                length.out = 1000)
)

# -----------------------------
# 3) Predictions + 95% CIs (approx. Bayesian)
# -----------------------------
pred_ci <- function(mod, xdat) {
  p <- predict(mod, newdata = xdat, type = "response", se.fit = TRUE)
  tibble(
    Measure = xdat$Measure,
    fit = p$fit,
    lo  = p$fit - 1.96 * p$se.fit,
    hi  = p$fit + 1.96 * p$se.fit
  )
}

pred_ba <- pred_ci(m_ba, newx) %>% mutate(elem = "Ba/Ca")
pred_mg <- pred_ci(m_mg, newx) %>% mutate(elem = "Mg/Ca")
pred_sr <- pred_ci(m_sr, newx) %>% mutate(elem = "Sr/Ca")

pred_all <- bind_rows(pred_ba, pred_mg, pred_sr)

# -----------------------------
# 4) Scaling and secondary axis
#    (Sr/Ca rescaled to fit left axis)
# -----------------------------
left_max_raw  <- max(data$Ba138, data$Mg24, pred_ba$hi, pred_mg$hi, na.rm = TRUE)
right_max_raw <- max(data$Sr87,  pred_sr$hi, na.rm = TRUE)
scale_sr <- left_max_raw / right_max_raw

# Long format + Sr/Ca rescaling for plotting only
data_long <- data %>%
  pivot_longer(cols = c(Ba138, Mg24, Sr87),
               names_to = "elem_raw", values_to = "value") %>%
  mutate(
    elem = recode(elem_raw,
                  "Ba138" = "Ba/Ca",
                  "Mg24"  = "Mg/Ca",
                  "Sr87"  = "Sr/Ca"),
    value_plot = ifelse(elem == "Sr/Ca", value * scale_sr, value)
  )

pred_all <- pred_all %>%
  mutate(
    fit_plot = ifelse(elem == "Sr/Ca", fit * scale_sr, fit),
    lo_plot  = ifelse(elem == "Sr/Ca", lo  * scale_sr, lo),
    hi_plot  = ifelse(elem == "Sr/Ca", hi  * scale_sr, hi)
  )

# Pretty breaks for each axis
left_breaks  <- pretty(c(0, left_max_raw),  n = 6)
left_limit   <- range(left_breaks)
right_breaks <- pretty(c(0, right_max_raw), n = 6)

# -----------------------------
# 5) Vertical lines (annual rings)
# -----------------------------
vlines <- c(0, 633.4, 801.4, 930.7, 1072.9, 1202.1, 1292.6,
            1409.0, 1499.4, 1615.8, 1758.0)

years_labels <- c(0, 1, 2, 3, 4,
                  5, 6, 7, 8, 9, 10)

# Shift labels slightly to the right of the lines
labels_x <- vlines + 40

# -----------------------------
# 6) Color palette
# -----------------------------
pal <- c("Ba/Ca" = "#E41A1C",  # red
         "Mg/Ca" = "#4DAF4A",  # green
         "Sr/Ca" = "#40BBD2")  # cyan

# -----------------------------
# 7) Plot
# -----------------------------
y_top <- max(left_limit, na.rm = TRUE) * 1.02

ggplot() +
  # Observed points (NO legend)
  geom_point(
    data = data_long,
    aes(Measure, value_plot, color = elem),
    alpha = 0.35, size = 1.4, show.legend = FALSE
  ) +
  
  # 95% ribbons (NO legend)
  geom_ribbon(
    data = pred_all,
    aes(Measure, ymin = lo_plot, ymax = hi_plot, fill = elem),
    alpha = 0.20, show.legend = FALSE
  ) +
  
  # Fitted curves (WITH legend)
  geom_line(
    data = pred_all,
    aes(Measure, fit_plot, color = elem),
    linewidth = 1.1, show.legend = TRUE
  ) +
  
  # Vertical lines (annual rings)
  geom_vline(xintercept = vlines, linetype = "dashed", linewidth = 0.4, alpha = 0.5) +
  
  # "Core" | numbers | "Edge" labels
  annotate("text",
           x = min(data$Measure, na.rm = TRUE), y = y_top,
           label = "Core", hjust = 0, vjust = 0, size = 4.2, fontface = "bold") +
  annotate("text",
           x = max(data$Measure, na.rm = TRUE), y = y_top,
           label = "Edge", hjust = 1, vjust = 0, size = 4.2, fontface = "bold") +
  annotate("text", x = labels_x, y = left_max_raw*1.02,
           label = 0:10, size = 3.5) +

  
  # Axes
  scale_y_continuous(
    name = expression("Ba:Ca, Mg:Ca  (mg g"^{-1}*")"),
    breaks = left_breaks,
    limits = left_limit,
    sec.axis = sec_axis(
      ~ . / scale_sr,
      name   = expression("Sr:Ca  (mg g"^{-1}*")"),
      breaks = right_breaks
    )
  ) +
  scale_x_continuous(
    name = expression("Core to otolith edge  ("*mu*"m)")
  ) +
  
  scale_color_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",  # legend outside panel
    legend.justification = "top",
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
    legend.title = element_blank(),
    plot.margin = margin(10, 25, 10, 10),
    axis.title.y.right = element_text(margin = margin(l = 10))
  )





###############################
# END
###############################
# Patrick Derviche
