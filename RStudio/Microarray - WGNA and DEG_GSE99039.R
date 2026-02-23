# ==============================================================================
# Script: PD Analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Package Installation & Loading
# ------------------------------------------------------------------------------

# Load libraries
library(GEOquery)
library(oligo)
library(Biobase)
library(splitstackshape)
library(tidyr)
library(dplyr)
library(arrayQualityMetrics)
library(tidyverse)
library(affy)
library(ggplot2)
library(CorLevelPlot)
library(gridExtra)
library(pROC)
library(limma)
library(tibble)
library(vegan)
library(ggpubr)
library(rstatix)
library(WGCNA)

# WGCNA Options: Important configuration.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA.
enableWGCNAThreads()

# ------------------------------------------------------------------------------
# 2. Data Loading & Metadata Preparation
# ------------------------------------------------------------------------------

# Exploring metadata
GSE99039 <- getGEO('GSE99039', GSEMatrix = TRUE)
length(GSE99039)
class(GSE99039[[1]])

# Prepare and analyze metadata labels:
varLabels(GSE99039[[1]]) 
info <- pData(GSE99039[[1]]) 

# Information
head(info)
# write.csv(info, file = "info_metadata.csv")

# Select columns of interest (1, 2, 50, 51, 57)
info.mod <- info[,c(1,2,50,51,57)]
colnames(info.mod)

# Modify the name to keep only the sample ID
names(info.mod) <- gsub(':ch1', '', names(info.mod))
# write.csv(info.mod, file = "info.mod.csv")

# Analyze disease levels: there are 14 types. We are interested in: "CONTROL", "IPD"
levels(as.factor(info.mod$`disease label`))

# Save the experiment in an object 
GSE99039_obj <- GSE99039[[1]]

# ------------------------------------------------------------------------------
# 3. Raw Data Processing & Normalization
# ------------------------------------------------------------------------------

# Get supplementary files (RAW data)
getGEOSuppFiles("GSE99039")

# Untar files
untar("GSE99039/GSE99039_RAW.tar", exdir = 'data/')

# Reading in .cel files and attach metadata
gse <- ReadAffy(celfile.path = "data/", phenoData=phenoData(GSE99039_obj))

# Preliminary data review for DEG and WGCNA
summary(exprs(gse))[,1:200] 
boxplot(exprs(gse)[,1:200], outline=FALSE, ylab = "Expression level", xlab = "Samples", col = "red")
abline(h=median(exprs(gse)),col="blue")

# Check data scale. Usually logarithmic (0 to 16).
plotDensities(exprs(gse), legend = FALSE, main = "Densities of raw expression")
hist(exprs(gse), col="blue", border="white", breaks=100, main = "Histogram of raw expression")

# Data is not normalized. Perform RMA normalization.
normalized.data <- affy::rma(gse)

# Check if normalized
boxplot(normalized.data[1:200], outline=FALSE, ylab = "Expression level", xlab = "Samples", col = "red")
abline(h=median(exprs(normalized.data)),col="blue")
plotDensities(exprs(normalized.data), legend = FALSE, main = "Densities of normalized expression")
hist(exprs(normalized.data), col="blue", border="white", breaks=100, main = "Histogram of normalized expression")

# Create dataframe with normalized data for WGCNA
normalized.expr <- exprs(normalized.data) 
normalized.expr.dataframe <- as.data.frame(exprs(normalized.data))

# Save the original dataframe
# write.csv(normalized.expr.dataframe, file = "normalized.exprs.dataframe.csv") 

# ------------------------------------------------------------------------------
# 4. Annotation & Cleaning
# ------------------------------------------------------------------------------

# Map probe IDs with gene symbol
feature.data <- GSE99039$GSE99039_series_matrix.txt.gz@featureData@data
feature.data.selection <- feature.data[,c(1,2,11)]

# Modify data
colnames(feature.data.selection)[3] <- "SYMBOL"
feature.data.selection <- feature.data.selection %>% 
  mutate(SYMBOL = gsub("///.*", "", SYMBOL))

normalized.expr.dataframe <- normalized.expr.dataframe %>%
  rownames_to_column(var = "ID") %>%
  inner_join(.,feature.data.selection, by = "ID")
colnames(normalized.expr.dataframe)

# Save as csv file
# write.csv(normalized.expr.dataframe, "normalized.expr.annotated.csv")

# Modify dataframe to have SYMBOL as ID and remove probes
data <- normalized.expr.dataframe %>% 
  select(-c(ID, GB_ACC)) %>% 
  relocate("SYMBOL")

data$SYMBOL[data$SYMBOL == ""] <- NA
data <- na.omit(data)
colSums(is.na(data))

# Modify and remove duplicate genes 
head(data)[1:10]
has_rownames(data)
data.mod <- remove_rownames(data) 
has_rownames(data.mod)

# Dplyr does not work with duplicated data. Remove them.
duplicado <- duplicated(data.mod$SYMBOL)
sum(duplicado)
data.mod.uniq <- data.mod %>% 
  distinct(SYMBOL, .keep_all = TRUE)
sum(duplicated(data.mod.uniq$SYMBOL))

# Modify titles to keep only GSM and set SYMBOL as rownames
data.mod.uniq[1:10,1:10]
data.mod.uniq <- data.mod.uniq %>% 
  gather(key = 'samples', value = 'counts', -SYMBOL) %>% 
  mutate(samples = gsub("_.*", "", samples)) %>% 
  spread(key = 'samples', value = 'counts') %>% 
  column_to_rownames(var = 'SYMBOL')

# Save
# write.csv(data.mod.uniq, "data.mod.uniq.WGCNA.csv")

# ------------------------------------------------------------------------------
# 5. Filtering Samples (IPD vs Control)
# ------------------------------------------------------------------------------

# In this data, there are samples other than IPD and Control. Remove them.

data.mod.t <- t(data.mod.uniq)

ncol(data.mod.uniq)
nrow(data.mod.uniq)
ncol(data.mod.t)
nrow(data.mod.t)

data.mod.t.prefil <- as.data.frame(cbind(data.mod.t, info.mod[4])) # Add info.mod columns
data.mod.t.prefil <- data.mod.t.prefil %>% relocate("disease label")

ncol(data.mod.t.prefil)
nrow(data.mod.t.prefil)

data.mod.t.filter <- data.mod.t.prefil %>% 
  filter(data.mod.t.prefil$`disease label` == 'CONTROL' | data.mod.t.prefil$`disease label` == 'IPD')

data.mod.t.final <- select(data.mod.t.filter, -"disease label")

ncol(data.mod.t.final)
nrow(data.mod.t.final) 

# ------------------------------------------------------------------------------
# 6. Quality Control (WGCNA/DEG)
# ------------------------------------------------------------------------------

# Outlier analysis. In WGCNA rows must be samples and columns genes.
# Detect gene outliers.
gsg <- goodSamplesGenes(data.mod.t.final)
summary(gsg)
gsg$allOK
# If TRUE, all data passed. If FALSE, there are outliers.

# Detect sample outliers - hierarchical clustering - method 1
htree <- hclust(dist(data.mod.t.final), method = "average")
plot(htree, hang = -1, cex = 0.6)  

plot(htree, hang = -1, cex = 0.6)
rect.hclust(htree , k = 4, border = 2:6)
abline(h = 4, col = 'red')

# PCA - method 2
pca <- prcomp(data.mod.t.final)
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() + xlim(-100,100) + ylim(-60,60) +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# plotMDS - method 3 - DEG

# Disease types:
levels(factor(normalized.data$`disease label:ch1`))

# Plot MDS for all samples (checking for batch effects or grouping)
plotMDS(exprs(normalized.data), 
        col = c("blue", "red", "green", "yellow","orange", "grey", "black", "pink","aquamarine","tan","khaki","forestgreen","blueviolet","bisque1" )[as.factor(normalized.data$`disease label:ch1` == "CONTROL")], pch = 19)
legend("topright", 
       pch = 19, 
       col = c("blue", "red", "green", "yellow","orange", "grey", "black", "pink","aquamarine","tan","khaki","forestgreen","blueviolet","bisque"), 
       cex = 1,
       legend = levels(as.factor(normalized.data$`disease label:ch1`)))

plotMDS(exprs(normalized.data), 
        col = c("blue", "red")[as.factor(normalized.data$`disease label:ch1` == "CONTROL" & normalized.data$`disease label:ch1` == "IPD")], pch = 19)

dists <- vegdist(normalized.data, method = "bray", k = 2)
nmds <- metaMDS(dists)
stressplot(nmds)
str(nmds)

# Note: In this case, we have a high number of samples and they look very similar according to PCA. 
# Nothing is removed.

# ------------------------------------------------------------------------------
# 7. WGCNA Analysis
# ------------------------------------------------------------------------------

# WGCNA analysis requires normalized data. (Data is already normalized).
data.mod.t.final # already transposed

# --- Network Construction ---

# Choose a threshold
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(data.mod.t.final,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

# Visualization to pick power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.9, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# Choice: Power 12 chosen as definitive.

# Prepare numeric data for WGCNA
data.mod.t.final[] <- sapply(data.mod.t.final, as.numeric)

soft_power <- 14
temp_cor <- cor
cor <- WGCNA::cor

# Memory estimate w.r.t blocksize
bwnet <- blockwiseModules(data.mod.t.final,
                          maxBlockSize = 24000, # Appropriate for 4-8GB RAM.
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor

# --- Module Eigengenes ---
# Save information
module_eigengenes <- bwnet$MEs
# Print out a preview
head(module_eigengenes)

# Get number of genes for each module
table(bwnet$colors)

module <- table(bwnet$colors)
write.table(module, "Modulos_Soft_14.txt", sep = "\t", row.names = FALSE)

# Save more variables. Another way to obtain the above.
moduleLabels <- bwnet$colors
moduleColors <- labels2colors(bwnet$colors)

# See module labels
bwLabels = matchLabels(bwnet$colors, moduleLabels, pThreshold = 1e-7);
bwColors = labels2colors(bwLabels)
table(bwLabels)

# Number of blocks (genes)
nBlocks = length(bwnet$dendrograms)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

table(bwnet$colors)
table(bwnet$unmergedColors)

# Grey module = all genes that doesn't fall into other modules were assigned to the grey module

# ------------------------------------------------------------------------------
# 8. Module-Trait Relationship
# ------------------------------------------------------------------------------

# Transform categorical variables (e.g., Control, AD, MCI) into binary variables.
# Modify data to remove extra samples (matching rows)

traits <- info.mod %>% 
  filter(data.mod.t.prefil$`disease label` == 'CONTROL' | data.mod.t.prefil$`disease label` == 'IPD') %>%
  mutate(PD = ifelse(grepl("IPD", `disease label`), 1, 0)) %>% 
  mutate(CONTROL = ifelse(grepl("CONTROL", `disease label`), 1, 0)) %>% 
  mutate(Sex_bin = ifelse(grepl('Female', Sex), 1, 0)) %>% 
  select(6,7)

# Define numbers of genes and samples
nSamples <- nrow(data.mod.t.final)
nGenes <- ncol(data.mod.t.final)

# Calculate correlation between module genes and traits
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
# Calculate p-values for these correlations
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# Visualize module-trait association as a heatmap
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

# Cannot have categorical columns, convert row.names column to row names
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')
str(heatmap.data)
colnames(heatmap.data)

# In x put all column names from Sex_bin. For y put names of genes/modules
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[10:11],
             titleX = "PD",
             y = names(heatmap.data)[1:8],
             titleY = "Module",
             col = c("blue1", "skyblue", "white" , "pink", "red"),
             main = "Module-Disease Relationships") # According to x

plotEigengeneNetworks(heatmap.data, "",
                      plotDendrograms = FALSE, plotHeatmaps = TRUE,
                      marDendro = c(1,6,1,6), 
                      marHeatmap = c(7,7,1,2), 
                      cex.lab = 0.9, xLabelsAngle = 90)

# Extract genes from interesting modules as a list.
module.gene.mapping <- as.data.frame(bwnet$colors)

blue_list <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'blue') %>% 
  rownames()
brown_list <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'brown') %>% 
  rownames()

write.table(blue_list, "Module_blue_Soft_14.txt", sep = "\t", row.names = FALSE)
write.table(brown_list, "Module_brown_Soft_12.txt", sep = "\t", row.names = FALSE)

# ------------------------------------------------------------------------------
# 8B. Intramodular Analysis: Identify Highly Connected Modular Genes
# ------------------------------------------------------------------------------

# Calculate the module membership and the associated p-values
# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 

module.membership.measure <- cor(module_eigengenes, data.mod.t.final, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:10,1:10]

# Calculate the gene significance and associated p-values
gene.signf.corr <- cor(data.mod.t.final, traits$diagnosis_bin, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)

# ------------------------------------------------------------------------------
# 9. Differential Expression Analysis (DEG)
# ------------------------------------------------------------------------------

# Specify experimental design. Design experimental matrix.
# Start directly from Large ExpressionSet object

experimental.design <- model.matrix( ~ 0 + normalized.data[['disease label:ch1']])
colnames(experimental.design) <- levels(as.factor(normalized.data[['disease label:ch1']]))
# Ensure names match your data structure, checking for 14 types as commented
colnames(experimental.design) <- c("APD", "CBD", "CONTROL", "DRD", "DRDDYT", "GU", "GPD", "HD", "HDHDB", "IPD", "MSA", "PDDEM", "PSP", "VD")

contrast_matrix <- makeContrasts(IPD - CONTROL, levels=experimental.design)
contrast_matrix

# Find difference between disease group and control 
# Compensate outliers with empirical weights
aw <- arrayWeights(exprs(normalized.data), experimental.design)

linear.fit <- lmFit(normalized.data,experimental.design, weights = aw)
contrast.linear.fit <- contrasts.fit(linear.fit,contrasts=contrast_matrix)
contrast.results <- eBayes(contrast.linear.fit)
summary(decideTests(contrast.results,lfc=0.2))

# logFC: log2 fold change
# AveExpr: Average expression
# t: logFC divided by SE
# P.Value: Raw p-value
# adj.P.Val: Adjusted p-value

topTable(contrast.results)

nrow(normalized.data)

# Extract data from condition 
Condition_Final <- topTable(contrast.results, number = 54675, coef = 1, sort.by = "logFC" )
head(Condition_Final)

write.csv(Condition_Final, "Condition_Final_Raw.csv")

# Annotate gene IDs
Annotated_Condition_Final <- Condition_Final %>%
  rownames_to_column(var = "ID") %>% 
  inner_join(.,feature.data.selection, by = "ID")

Annotated_Condition_Final <- select(Annotated_Condition_Final, -GB_ACC)
Annotated_Condition_Final$SYMBOL[Annotated_Condition_Final$SYMBOL == ""] <- NA 

write.csv(Annotated_Condition_Final, "Annotated_Condition_Final_NA.csv")

colSums(is.na(Annotated_Condition_Final))
Annotated_Condition_Final <- na.omit(Annotated_Condition_Final)
colSums(is.na(Annotated_Condition_Final))

# Remove duplicate genes
head(Annotated_Condition_Final)
has_rownames(Annotated_Condition_Final)
Annotated_Condition_Final.mod <- remove_rownames(Annotated_Condition_Final) 
has_rownames(Annotated_Condition_Final.mod)

duplicado <- duplicated(Annotated_Condition_Final.mod$SYMBOL)
sum(duplicado)

Annotated_Condition_Final.mod <- Annotated_Condition_Final.mod %>% 
  distinct(SYMBOL, .keep_all = TRUE)
sum(duplicated(Annotated_Condition_Final.mod$SYMBOL))

# Save
write.csv(Annotated_Condition_Final, "Annotated_Condition_Final.csv")
write.csv(Annotated_Condition_Final.mod, "Annotated_Condition_Final.mod.csv")

# ------------------------------------------------------------------------------
# 10. Volcano Plot Visualization
# ------------------------------------------------------------------------------

# To work with this package, must work from a dataframe
class(Annotated_Condition_Final.mod)

volplot <- ggplot(Annotated_Condition_Final.mod, aes(x = logFC, y= -log10(adj.P.Val))) + geom_point() 

volplot +  geom_hline(yintercept = -log10(0.05),linetype = "dashed") + 
  geom_vline(xintercept = c(log2(1.07), log2(0.93)),linetype = "dashed") + 
  xlim(-1.5, 1.5)

# Create new categorical column
Annotated_Condition_Final_DEG <- Annotated_Condition_Final.mod %>%
  mutate(gene_type = case_when(Annotated_Condition_Final.mod$logFC > 0.1 & Annotated_Condition_Final.mod$adj.P.Val < 0.05 ~ "Up",
                               Annotated_Condition_Final.mod$logFC < -0.1 & Annotated_Condition_Final.mod$adj.P.Val < 0.05 ~ "Down",
                               TRUE ~ "NS"))   

# Obtain gene_type counts            
Annotated_Condition_Final_DEG %>%
  count(gene_type)

# Save this data
write.csv(Annotated_Condition_Final_DEG, "Annotated_Condition_Final_DEG_0.1.csv")

# Representation of genes in the plot
Annotated_Condition_Final_DEG %>%
  distinct(gene_type) %>%
  pull() 

# Add colour, size and alpha
cols <- c("Up" = "red", "Down" = "blue", "NS" = "grey") 
sizes <- c("Up" = 4, "Down" = 4, "NS" = 2) 
alphas <- c("Up" = 1, "Down" = 1, "NS" = 0.5)

mi_volcano_plot <- Annotated_Condition_Final_DEG %>%
  ggplot(aes(x = logFC, 
             y= -log10(adj.P.Val),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21,   
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(1.07), log2(0.93)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + xlim(-1.5, 1.5) + ylim(0,7) + # Modify point transparency and scale of axis 
  labs(title = NULL,
       x = expression(log[2]("Fold Change")),
       y = expression(-log[10]("Adjusted P-value")),
       fill = "Genes:",
       size = "Genes:",
       alpha = "Genes:") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        legend.position = "right") 

mi_volcano_plot

ggsave(
  filename = "Volcano_Plot_Final.tiff", 
  plot = mi_volcano_plot,               
  width = 12,                             
  height = 9,                           
  dpi = 600                             
)

# ------------------------------------------------------------------------------
# 11. Expression Analysis: T-Test
# ------------------------------------------------------------------------------

ncol(data.mod.t.filter)
nrow(data.mod.t.filter)

data.ttest <- data.mod.t.filter %>% 
  rename(diagnosis = `disease label`) %>% 
  filter(diagnosis == "CONTROL" | diagnosis == "IPD") %>% 
  select(diagnosis, LILRB2,TLN1,ITGAM,C5AR1,IL18,PTPRC,PTEN,CTSS,TNFSF13B,TLR4,CSF1R,NCF4,MYO1F,VAV1,HCLS1)
write.csv(data.ttest, "data.ttest_PD_original.csv", row.names = FALSE)

data.ttest %>% 
  group_by(diagnosis) %>% 
  get_summary_stats(HCLS1, type = "mean_sd")
leveneTest(HCLS1 ~ diagnosis, data = data.ttest) 
# If p value < 0.05 = variances are significantly different; If p value > 0.05 = homogeneous variances

data.ttest %>%
  group_by(diagnosis) %>%
  identify_outliers(HCLS1) # Check for outliers

df_HCLS1 <- data.ttest %>%
  group_by(diagnosis) %>% 
  mutate(Q1 = quantile(HCLS1, 0.25),      # First quartile (Q1)
         Q3 = quantile(HCLS1, 0.75),      # Third quartile (Q3)
         IQR = Q3 - Q1,                   # Calculate IQR
         Lower_Bound = Q1 - 1.5 * IQR,    # Lower bound
         Upper_Bound = Q3 + 1.5 * IQR) %>%
  filter(HCLS1 >= Lower_Bound & HCLS1 <= Upper_Bound) %>%  # Filter out outliers
  select(-Q1, -Q3, -IQR, -Lower_Bound, -Upper_Bound)       # Remove auxiliary columns

df_HCLS1 <- as.data.frame(df_HCLS1)
df_HCLS1 %>% 
  group_by(diagnosis) %>% 
  get_summary_stats(HCLS1, type = "mean_sd")
leveneTest(HCLS1 ~ diagnosis, data = df_HCLS1)

bxp_HCLS1 <- ggboxplot(
  df_HCLS1, x = "diagnosis", y = "HCLS1", 
  ylab = "HCLS1", xlab = "Diagnosis", add = "jitter")
bxp_HCLS1

stat.test_HCLS1 <- df_HCLS1 %>% 
  t_test(HCLS1 ~ diagnosis, var.equal = TRUE) %>%
  add_significance()
stat.test_HCLS1

stat.test_HCLS1 <- stat.test_HCLS1 %>% add_xy_position(x = "diagnosis")
x_HCLS1  <- bxp_HCLS1 + 
  stat_pvalue_manual(stat.test_HCLS1, tip.length = 0) +
  labs(subtitle = get_test_label(stat.test_HCLS1, detailed = TRUE))
x_HCLS1 
ggsave("ttest_HCLS1.png", plot = x_HCLS1, width = 12, height = 10, dpi = 300)