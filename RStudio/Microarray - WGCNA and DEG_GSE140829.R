# ==============================================================================
# Script: AD Analysis
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Package Installation & Loading
# ------------------------------------------------------------------------------

# Install packages (Run only if necessary)
# BiocManager::install("affy")
# BiocManager::install("oligo")
# BiocManager::install("Biobase")
# BiocManager::install("GEOquery")
# BiocManager::install("arrayQualityMetrics")
# install.packages("remotes")
# remotes::install_github("kevinblighe/CorLevelPlot")

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
library(ggpubr)
library(rstatix)
library(car)
library(WGCNA)

# WGCNA Options: Important configuration.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA (skip if using RStudio/third-party environments if it causes issues).
enableWGCNAThreads()

# ------------------------------------------------------------------------------
# 2. Data Loading & Metadata Preparation
# ------------------------------------------------------------------------------


# Exploring metadata
GSE140829 <- getGEO('GSE140829', GSEMatrix = TRUE)
length(GSE140829)
class(GSE140829[[1]])

# Prepare and analyze metadata labels:
varLabels(GSE140829[[1]]) # View labels associated with this experiment
info <- pData(GSE140829[[1]]) ## Print the sample information

# Information
head(info)
# write.csv(info, file = "info_metadata.csv")

# Select columns of interest (1, 2, 41, 42)
info <- info[,c(1,2,41:42)]
colnames(info)

# Modify the name to keep only the sample ID
info.mod <- 
  info %>%
  mutate(title = gsub("Whole blood,", " ", title)) %>%
  mutate(title = gsub(".*,", " ", title)) %>%
  mutate(title = gsub("\\[|\\]", " ", title)) %>%
  mutate(title = gsub("ad_mci", " ", title)) %>%
  as.data.frame()
names(info.mod) <- gsub(':ch1', '', names(info.mod))
# write.csv(info.mod, file = "info.mod.csv")

# ------------------------------------------------------------------------------
# 3. Data Preprocessing & Annotation
# ------------------------------------------------------------------------------

# Save expression data into an object for DEG analysis. 
# Starting from .txt data. Save data into gse object.
# Object Large ExpressionSet
gse <- GSE140829[[1]]

# Preliminary data review for DEG and WGCNA
summary(exprs(gse))[,1:200] # Check data expression
boxplot(exprs(gse)[,1:200], outline=FALSE, ylab = "Expression level", xlab = "Samples", col = "red")
abline(h=median(exprs(gse)),col="blue")
# This allows checking the data scale. Usually found in a logarithmic scale (0 to 16).

# Create a dataframe from this object for saving and WGCNA
gse.expr <- exprs(gse)
normalized.expr.dataframe <- as.data.frame(gse.expr)

# Save the original dataframe
# write.csv(normalized.expr.dataframe, file = "normalized.exprs.dataframe.csv") 

# Add each probe with its SYMBOL
feature.data <- GSE140829$GSE140829_series_matrix.txt.gz@featureData@data
feature.data.selection <- feature.data[,c(1,3)]
colnames(feature.data.selection)[2] <- "SYMBOL"

# Join dataframes
normalized.expr.dataframe.ann <- normalized.expr.dataframe %>%
  rownames_to_column(var = "ID") %>%
  inner_join(.,feature.data.selection, by = "ID") 

# Check
colnames(normalized.expr.dataframe.ann)
# write.csv(normalized.expr.dataframe.ann, "normalized.expr.annotated.csv")

# Modify dataframe to have SYMBOL as ID and remove probes
data <- normalized.expr.dataframe.ann %>% 
  select(-c(ID)) %>% 
  relocate("SYMBOL") %>% 
  na.omit() 

# Modify, remove gene duplicates and set SYMBOL as rowname
head(data)[1:10]
has_rownames(data)
data.mod <- remove_rownames(data) 
has_rownames(data.mod)

# Dplyr does not work with duplicated data. They must be removed.
duplicado <- duplicated(data.mod$SYMBOL)
sum(duplicado)
data.mod.uniq <- data.mod %>% 
  distinct(SYMBOL, .keep_all = TRUE)
sum(duplicated(data.mod.uniq$SYMBOL))
data.mod.uniq <- column_to_rownames(data.mod.uniq, var = 'SYMBOL')

# Save
# write.csv(data.mod, "data.mod.uniq.WGCNA.csv")

# ------------------------------------------------------------------------------
# 4. Filtering Samples (AD vs Control)
# ------------------------------------------------------------------------------

# In this data, there are not only AD and Control, but also MCI. We will remove them.

data.mod.t <- t(data.mod.uniq)

ncol(data.mod.uniq)
nrow(data.mod.uniq)
ncol(data.mod.t)
nrow(data.mod.t)

data.mod.t.filter <- as.data.frame(cbind(data.mod.t, info.mod[3])) # Add info.mod columns
data.mod.t.filter <- data.mod.t.filter %>% relocate("diagnosis")

ncol(data.mod.t.filter)
nrow(data.mod.t.filter)

# Filter to keep only Control and AD
data.mod.t.final <- data.mod.t.filter %>% 
  filter(data.mod.t.filter$diagnosis == 'Control' | data.mod.t.filter$diagnosis == 'AD')

data.mod.t.final <- select(data.mod.t.final, -diagnosis)

ncol(data.mod.t.final)
nrow(data.mod.t.final) 

# ------------------------------------------------------------------------------
# 5. Quality Control (WGCNA/DEG)
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
  geom_point() + xlim(-60,30) + ylim(-30,30) +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# plotMDS - method 3 - DEG
varLabels(gse)
pData(gse) # dataset information

# Our information of main interest is "diagnosis:ch1"
# Analyze outliers to see possible contamination and if groups are clearly differentiated:

plotMDS(exprs(gse), labels=gse$`diagnosis:ch1`, pch=20)
plotMDS(exprs(gse), 
        col = c("blue", "green", "red")[as.factor(gse$`diagnosis:ch1`)], pch = 19)
legend("topright", 
       pch = 19, 
       col = c("blue", "green", "red"), 
       cex = 1,
       legend = levels(as.factor(gse$`diagnosis:ch1`)))

# If batch effect is observed, it must be corrected.
# In this case, we have a high number of samples and they are very similar. Nothing is removed.

# ------------------------------------------------------------------------------
# 6. WGCNA Analysis
# ------------------------------------------------------------------------------

# WGCNA analysis requires normalized data (already normalized).
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
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# Choice: Power 8 chosen (based on comments).

# Prepare numeric data for WGCNA
data.mod.t.final[] <- sapply(data.mod.t.final, as.numeric)

soft_power <- 12 # Note: Comment said 8, code uses 12. Adjust as needed.
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
write.table(module, "Modulos_Soft_12.txt", sep = "\t", row.names = FALSE)

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

# Plot the dendrogram and the module colors underneath for each block
for (block in 1:nBlocks) {
  plotDendroAndColors(bwnet$dendrograms[[block]], moduleColors[bwnet$blockGenes[[block]]],
                      "Module colors",
                      main = paste("Gene dendrogram and module colors in block", block),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      setLayout = FALSE)
}

# Grey module = all genes that doesn't fall into other modules were assigned to the grey module

# ------------------------------------------------------------------------------
# 7. Module-Trait Relationship
# ------------------------------------------------------------------------------

# Transform categorical variables (e.g., Control, AD, MCI) into binary variables.
# Modify data to remove MCI from traits (to match rows)

traits <- info.mod %>% 
  filter(data.mod.t.filter$diagnosis == 'Control' | data.mod.t.filter$diagnosis == 'AD') %>% 
  mutate(AD = ifelse(grepl("AD", diagnosis), 1, 0)) %>% 
  mutate(CONTROL = ifelse(grepl("Control", diagnosis), 1, 0)) %>% 
  mutate(Sex_bin = ifelse(grepl('Female', Sex), 1, 0)) %>% 
  select(5,6)

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

# In x put all column names from Sex_bin. For y put names of genes/modules
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[21:22],
             titleX = "AD",
             y = names(heatmap.data)[1:19],
             titleY = "Module",
             col = c("blue1", "skyblue", "white" , "pink", "red"),
             main = "Module-Disease Relationships") # According to x

plotEigengeneNetworks(heatmap.data, "",
                      plotDendrograms = FALSE, plotHeatmaps = TRUE,
                      marDendro = c(1,6,1,6), 
                      marHeatmap = c(7,7,1,2), 
                      cex.lab = 0.9, xLabelsAngle = 90)

# Higher number of asterisks = trait is significantly associated with that module.
# Higher value (more positive = redder) = more correlated with value assigned to 1.

# Extract genes from interesting modules as a list.
module.gene.mapping <- as.data.frame(bwnet$colors)

red_list <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'red') %>% 
  rownames()
write.table(red_list, "Module_red_12.txt", sep = "\t", row.names = FALSE) # Fixed variable name

magenta_list <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'magenta') %>% 
  rownames()
write.table(magenta_list, "Module_magenta_12.txt", sep = "\t", row.names = FALSE)

pink_list <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'pink') %>% 
  rownames()
write.table(pink_list, "Module_pink_12.txt", sep = "\t", row.names = FALSE)

# ------------------------------------------------------------------------------
# 7B. Intramodular Analysis: Identify Highly Connected Modular Genes
# ------------------------------------------------------------------------------

# Calculate the module membership and the associated p-values
# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 

module.membership.measure <- cor(module_eigengenes, data.mod.t.final, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:10,1:10]

# Calculate the gene significance and associated p-values
gene.signf.corr <- cor(data.mod.t.final, traits$diagnosis_bin, use = 'p') # Check if diagnosis_bin exists in traits
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)

# ------------------------------------------------------------------------------
# 8. Differential Expression Analysis (DEG)
# ------------------------------------------------------------------------------

# Specify experimental design.
varLabels(gse)

experimental.design <- model.matrix( ~ 0 + gse[["diagnosis:ch1"]])
colnames(experimental.design) <- levels(as.factor(gse[['diagnosis:ch1']]))
colnames(experimental.design) <- c("AD", "Control", "MCI")
experimental.design

contrast_matrix <- makeContrasts(AD - Control, levels=experimental.design)
contrast_matrix

# Find difference between disease group and control 
# Compensate outliers with empirical weights
aw <- arrayWeights(exprs(gse), experimental.design)

linear.fit <- lmFit(gse,experimental.design, weights = aw)
contrast.linear.fit <- contrasts.fit(linear.fit,contrasts=contrast_matrix)
contrast.results <- eBayes(contrast.linear.fit)
summary(decideTests(contrast.results,lfc=0.1))

# logFC: log2 fold change
# AveExpr: Average expression
# t: logFC divided by SE
# P.Value: Raw p-value
# adj.P.Val: Adjusted p-value

topTable(contrast.results)
table(decideTests(contrast.results, lfc=0.1))

# Check total number of genes
nrow(gse)

# Extract data from condition 
Condition_Final <- topTable(contrast.results, number = 47008, coef = 1, sort.by = "logFC" )
head(Condition_Final)

write.csv(Condition_Final, "Condition_Final_Raw.csv")

# Annotate gene IDs (Already present in Gene Symbol)
Annotated_Condition_Final <- Condition_Final[3:10]

# Clean empty or null data in Gene Symbol
Annotated_Condition_Final$GENE.SYMBOL[Annotated_Condition_Final$GENE.SYMBOL == ""] <- NA 
write.csv(Annotated_Condition_Final, "Annotated_Condition_Final_NA.csv")

colSums(is.na(Annotated_Condition_Final))

# Remove NA values
Annotated_Condition_Final <- na.omit(Annotated_Condition_Final)
colSums(is.na(Annotated_Condition_Final))

# Remove duplicate genes
head(Annotated_Condition_Final)
has_rownames(Annotated_Condition_Final)
Annotated_Condition_Final.mod <- remove_rownames(Annotated_Condition_Final) 
has_rownames(Annotated_Condition_Final.mod)

# Dplyr does not work with duplicated data. Remove them.
duplicado <- duplicated(Annotated_Condition_Final.mod$GENE.SYMBOL)
sum(duplicado)

Annotated_Condition_Final.mod <- Annotated_Condition_Final.mod %>% 
  distinct(GENE.SYMBOL, .keep_all = TRUE)
sum(duplicated(Annotated_Condition_Final.mod$GENE.SYMBOL))

# Save
write.csv(Annotated_Condition_Final.mod, "Annotated_Condition_Final.mod.csv")

# ------------------------------------------------------------------------------
# 9. Volcano Plot Visualization
# ------------------------------------------------------------------------------

# To work with ggplot, must work from a dataframe
class(Annotated_Condition_Final.mod)

volplot <- ggplot(Annotated_Condition_Final.mod, aes(x = logFC, y= -log10(adj.P.Val))) + geom_point() 

volplot +  geom_hline(yintercept = -log10(0.05),linetype = "dashed") + 
  geom_vline(xintercept = c(log2(1.07), log2(0.93)),linetype = "dashed") + 
  xlim(-0.5, 0.5)

# Plot limit is between -0.5 and 0.5. Ideal cutoff is 0.1 / -0.1.

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
  scale_alpha_manual(values = alphas) + xlim(-0.5, 0.5) + ylim(0,8) + # Modify point transparency and scale of axis 
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
# 10. Expression Analysis: T-Test
# ------------------------------------------------------------------------------

# Use dataframe data.mod.filter. Filter selecting two groups: AD and Control. 
# Select specific genes (e.g., HCLS1)

data.ttest <- data.mod.t.filter %>% 
  filter(diagnosis == "Control" | diagnosis == "AD") %>% 
  select(diagnosis, LILRB2,TLN1,ITGAM,C5AR1,RPL15,FOS,EIF4A2,NKTR,TPT1,RBM25,RPL22,TAX1BP1,RPS6KB1,NCF4,MYO1F,VAV1,HCLS1)
write.csv(data.ttest, "data.ttest_AD_original.csv", row.names = FALSE)

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
  t_test(HCLS1 ~ diagnosis, var.equal = FALSE) %>%
  add_significance()
stat.test_HCLS1

stat.test_HCLS1 <- stat.test_HCLS1 %>% add_xy_position(x = "diagnosis")
x_HCLS1  <- bxp_HCLS1 + 
  stat_pvalue_manual(stat.test_HCLS1, tip.length = 0) +
  labs(subtitle = get_test_label(stat.test_HCLS1, detailed = TRUE))
x_HCLS1 
ggsave("ttest_HCLS1.png", plot = x_HCLS1, width = 12, height = 10, dpi = 300)
