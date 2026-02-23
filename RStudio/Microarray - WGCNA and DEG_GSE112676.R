# ------------------------------------------------------------------------------
# Script: ALS Analysis
# 1. Package Installation & Loading
# ------------------------------------------------------------------------------

# Install packages (Run only if necessary)
BiocManager::install("affy")
BiocManager::install("oligo")
BiocManager::install("Biobase")
BiocManager::install("GEOquery")
BiocManager::install("arrayQualityMetrics")

# Additional packages to consider
# install.packages("splitstackshape")
# install.packages("tidyverse") 

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

# WGCNA Options: Important configuration, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA.

# ------------------------------------------------------------------------------
# 2. Data Loading & Metadata Preparation
# ------------------------------------------------------------------------------

# Exploring metadata
GSE112676 <- getGEO('GSE112676')
length(GSE112676)
class(GSE112676[[1]])

# Prepare and analyze metadata labels:
varLabels(GSE112676[[1]]) # View labels associated with this experiment
info <- pData(GSE112676[[1]]) ## Print the sample information

# Information
head(info)
# write.csv(info, file = "info_metadata.csv")

# Select columns of interest
info <- info[,c(1,2,42,43,45)]
colnames(info)

# Modify the name to keep only the sample ID
info.mod <- 
  info %>%
  mutate(title = gsub("whole blood gene expression profile_", "", title)) %>%
  as.data.frame()
names(info.mod) <- gsub(':ch1', '', names(info.mod))
write.csv(info.mod, file = "info.mod.csv")

# ------------------------------------------------------------------------------
# 3. Data Preprocessing & Annotation
# ------------------------------------------------------------------------------

# Save expression data into an object for DEG analysis. 
# Starting from .txt data. Save data into gse object.
# Object Large ExpressionSet
gse <- GSE112676[[1]]

# Preliminary data review for DEG and WGCNA
summary(exprs(gse))[,1:200] # Check data expression
boxplot(exprs(gse)[,1:200], outline=FALSE, ylab = "Expression level", xlab = "Samples", col = "red")
abline(h=median(exprs(gse)),col="blue")
# This allows checking the data scale.
# Usually found in a logarithmic scale (0 to 16).

# Create a dataframe from this object for saving and WGCNA
gse.expr <- exprs(gse)
normalized.expr.dataframe <- as.data.frame(gse.expr)

# Save the original dataframe
write.csv(normalized.expr.dataframe, file = "normalized.exprs.dataframe.csv")

# Add each probe with its SYMBOL
feature.data <- GSE112676$GSE112676_series_matrix.txt.gz@featureData@data
feature.data.selection <- feature.data[,c(1,14)]
colnames(feature.data.selection)[2] <- "SYMBOL"

# Join dataframes
normalized.expr.dataframe.ann <- normalized.expr.dataframe %>%
  rownames_to_column(var = "ID") %>%
  inner_join(.,feature.data.selection, by = "ID") 

# Check
colnames(normalized.expr.dataframe.ann)

# Save as csv file
write.csv(normalized.expr.dataframe.ann, "normalized.expr.annotated.csv")

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
write.csv(data.mod, "data.mod.uniq.WGCNA.csv")

# In this data, there are only ALS and Control samples. No modification needed.
data.mod.t <- t(data.mod.uniq)

ncol(data.mod.t)
nrow(data.mod.t)

# ------------------------------------------------------------------------------
# 4. Quality Control (WGCNA/DEG)
# ------------------------------------------------------------------------------

# Outlier analysis. In WGCNA rows must be samples and columns genes.

# Detect gene outliers.
gsg <- goodSamplesGenes(data.mod.t)
summary(gsg)
gsg$allOK
# If TRUE, all data passed. If FALSE, there are outliers.

# Detect sample outliers - hierarchical clustering - method 1
htree <- hclust(dist(data.mod.t), method = "average")

plot(htree, hang = -1, cex = 0.6)
rect.hclust(htree , k = 4, border = 2:6)
abline(h = 4, col = 'red')

# PCA - method 2
pca <- prcomp(data.mod.t)
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() + xlim(-50,60) + ylim(-50,50) +
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
        col = c("blue", "red")[as.factor(gse$`diagnosis:ch1`)], pch = 19)
legend("topright", 
       pch = 19, 
       col = c("blue", "red"), 
       cex = 1,
       legend = levels(as.factor(gse$`diagnosis:ch1`)))

# Do not forget: if a batch effect is seen, it must be corrected.
# In this case, we have a high number of samples and they are very similar. Nothing is removed.
# Example code for exclusion if needed:
# samples.to.be.excluded <- c('GSM4615000', 'GSM4614993', 'GSM4614995')
# data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]

# ------------------------------------------------------------------------------
# 5. WGCNA Analysis
# ------------------------------------------------------------------------------

# WGCNA analysis requires normalized data. (Data is already normalized).
# We have metadata and normalized gene expression data.

data.mod.t # already transposed

# --- Network Construction ---

# Choose a threshold
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(data.mod.t,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

# Visualization to pick power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.85, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# When choosing the threshold, we look for R^2 above 0.8 but with a low mean.
# Usually around 6. Pick the smallest power with the lowest possible mean value.

# To execute WGCNA block function, expression data needs to be numeric
data.mod.t[] <- sapply(data.mod.t, as.numeric)

soft_power <- 10
temp_cor <- cor
cor <- WGCNA::cor

# Memory estimate w.r.t blocksize
bwnet <- blockwiseModules(data.mod.t,
                          maxBlockSize = 20000, # Appropriate for 4-8GB RAM. Try to keep as 1 block (= nGenes)
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 4586,
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
write.table(module, "Modulos_Soft_10.txt", sep = "\t", row.names = FALSE)

# Number of blocks (genes)
nBlocks = length(bwnet$dendrograms)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05,
                    main = paste("Gene dendrogram and module colors"))

# Grey module = all genes that doesn't fall into other modules were assigned to the grey module

# ------------------------------------------------------------------------------
# 6. Module-Trait Relationship
# ------------------------------------------------------------------------------

# To analyze association between modules and phenotypic traits, transform categorical 
# variables (e.g., Control, ALS) into binary variables // Control = 0, ALS = 1 

traits <- info.mod %>% 
  mutate(ALS = ifelse(grepl("ALS", diagnosis), 1, 0)) %>% 
  mutate(CONTROL = ifelse(grepl("CON", diagnosis), 1, 0)) %>% 
  mutate(Sex_bin = ifelse(grepl('Female', Sex), 1, 0)) %>% 
  select(6, 7)

# Define numbers of genes and samples
nSamples <- nrow(data.mod.t)
nGenes <- ncol(data.mod.t)

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
             x = names(heatmap.data)[13:14],
             titleX = " ",
             y = names(heatmap.data)[1:11],
             titleY = " ",
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

b_colors <- as.data.frame(bwnet$colors)
write.csv(b_colors, "b_colors.csv")

greenyellow_list <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'greenyellow') %>% 
  rownames()
write.table(greenyellow_list, "Module_greenyellow_Soft_10.txt", sep = "\t", row.names = FALSE)

magenta_list <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'magenta') %>% 
  rownames()
write.table(magenta_list, "Module_magenta_Soft_10.txt", sep = "\t", row.names = FALSE)

# ------------------------------------------------------------------------------
# 6B. Intramodular Analysis: Identify Highly Connected Modular Genes
# ------------------------------------------------------------------------------

# Calculate the module membership and the associated p-values
# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, data.mod.t, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:3,1:3]

# Calculate the gene significance and associated p-values
gene.signf.corr <- cor(data.mod.t, traits$diagnosis_bin, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)

# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.

# ------------------------------------------------------------------------------
# 7. Differential Expression Analysis (DEG)
# ------------------------------------------------------------------------------

varLabels(gse) 

# Our information of main interest is "diagnosis:ch1"
ALS <- gse$`diagnosis:ch1` == "ALS"
CONTROL <- gse$`diagnosis:ch1` == "CON"
sum(ALS) 
sum(CONTROL) 

# Specify experimental design. Design experimental matrix.
# Contrast matrix. Must start from ExpressionSet object.

experimental.design <- model.matrix( ~ 0 + gse[['diagnosis:ch1']])
colnames(experimental.design) <- levels(as.factor(gse[['diagnosis:ch1']]))
colnames(experimental.design) <- c("ALS", "CONTROL") 
experimental.design

contrast_matrix <- makeContrasts(ALS - CONTROL, levels=experimental.design)
contrast_matrix

# Find difference between disease group and control 
# Created contrast studies. -1 is control, 0 condition does not act, 1 is condition to compare.

# Compensate outliers with empirical weights
aw <- arrayWeights(exprs(gse), experimental.design)

linear.fit <- lmFit(gse,experimental.design, weights = aw)
contrast.linear.fit <- contrasts.fit(linear.fit,contrasts=contrast_matrix)
contrast.results <- eBayes(contrast.linear.fit)
summary(decideTests(contrast.results,lfc=0.2))

# logFC: log2 fold change
# AveExpr: Average expression across all samples, in log2 CPM
# t: logFC divided by its standard error
# P.Value: Raw p-value
# adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
# B: log-odds that gene is DE

topTable(contrast.results)
# table(decideTests(contrast.results))

# To know total number of genes, check initial visualization or use function:
nrow(gse)

# Extract data from condition 
Condition_Final <- topTable(contrast.results, number = 47008, coef = 1, sort.by = "logFC" )
head(Condition_Final)

write.csv(Condition_Final, "Condition_Final_Raw.csv")

# Annotate gene IDs
# In this case, not applied because Gene Symbol IDs are already included.

Annotated_Condition_Final <- Condition_Final[,c(14,31,32,33,34,35,36)]

# Clean empty or null data in Gene Symbol and Entrez_Gene_ID columns
Annotated_Condition_Final$Symbol[Annotated_Condition_Final$Symbol == " "] <- NA

write.csv(Annotated_Condition_Final, "Annotated_Condition_Final_NA.csv")

# Changed empty values to NA. Analyze number of NA values per column.
colSums(is.na(Annotated_Condition_Final))

# Remove duplicate genes
# Modify, remove gene duplicates and set SYMBOL as rowname
head(Annotated_Condition_Final)
has_rownames(Annotated_Condition_Final)
Annotated_Condition_Final.mod <- remove_rownames(Annotated_Condition_Final) 
has_rownames(Annotated_Condition_Final.mod)

# Dplyr does not work with duplicated data. Remove them.
duplicado <- duplicated(Annotated_Condition_Final.mod$Symbol)
sum(duplicado)

Annotated_Condition_Final.mod <- Annotated_Condition_Final.mod %>% 
  distinct(Symbol, .keep_all = TRUE)
sum(duplicated(Annotated_Condition_Final.mod$GENE.SYMBOL))

# Save files:
write.csv(Annotated_Condition_Final, "Annotated_Condition_Final.csv")

# ------------------------------------------------------------------------------
# 8. Volcano Plot Visualization
# ------------------------------------------------------------------------------

# To work with this package, must work from a dataframe
class(Annotated_Condition_Final.mod)

volplot <- ggplot(Annotated_Condition_Final.mod, aes(x = logFC, y= -log10(adj.P.Val))) + geom_point()

volplot +  geom_hline(yintercept = -log10(0.05),linetype = "dashed") + 
  geom_vline(xintercept = c(log2(1.148), log2(0.87)),linetype = "dashed") + 
  xlim(-1, 1)

# Keep 0.2 cutoff for now

# Create new categorical column
Annotated_Condition_Final_DEG <- Annotated_Condition_Final.mod %>%
  mutate(gene_type = case_when(Annotated_Condition_Final.mod$logFC > 0.2 & Annotated_Condition_Final.mod$adj.P.Val < 0.05 ~ "Up",
                               Annotated_Condition_Final.mod$logFC < -0.2 & Annotated_Condition_Final.mod$adj.P.Val < 0.05 ~ "Down",
                               TRUE ~ "NS"))   
# Remember logFC gives values in log base 2, so threshold is 0 = log2(1)

# Obtain gene_type counts            
Annotated_Condition_Final_DEG %>%
  count(gene_type)

# Save this data
write.csv(Annotated_Condition_Final_DEG, "ALS_DEG.csv")

# Representation of genes in the plot
Annotated_Condition_Final_DEG %>%
  distinct(gene_type) %>%
  pull() 

# Add colour, size and alpha (transparency) to volcano plot
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
  geom_vline(xintercept = c(log2(1.148), log2(0.87)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + xlim(-1, 1) + ylim(0,60) + # Modify point transparency and scale of axis 
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
# 9. Expression Analysis: T-Test
# ------------------------------------------------------------------------------

# Use dataframe data.mod.filter. Filter selecting two groups: ALS and Control. Select genes
# U2AF2 and CHD4 (example genes mentioned)

data.mod.t.filter <- as.data.frame(cbind(data.mod.t, info.mod[3])) # Add info.mod columns
data.mod.t.filter <- data.mod.t.filter %>% relocate("diagnosis")

ncol(data.mod.t.filter)
nrow(data.mod.t.filter)

data.ttest <- data.mod.t.filter %>% 
  filter(diagnosis == "CON" | diagnosis == "ALS") %>% 
  select(diagnosis, LILRB2,TLN1,ITGAM,C5AR1,RPL15,FOS,EIF4A2,NKTR,TPT1,RBM25,RPL22,TAX1BP1,RPS6KB1,IL18,PTPRC,PTEN,CTSS,TNFSF13B,TLR4,CSF1R)
write.csv(data.ttest, "data.ttest_ALS_original.csv", row.names = FALSE)

data.ttest %>% 
  group_by(diagnosis) %>% 
  get_summary_stats(CSF1R, type = "mean_sd")
leveneTest(CSF1R ~ diagnosis, data = data.ttest) 
# If p value < 0.05 = variances are significantly different; If p value > 0.05 = homogeneous variances
# If homogeneous, ok for t-test

data.ttest %>%
  group_by(diagnosis) %>%
  identify_outliers(CSF1R) 

df_CSF1R <- data.ttest %>%
  group_by(diagnosis) %>% 
  mutate(Q1 = quantile(CSF1R, 0.25),      # First quartile (Q1)
         Q3 = quantile(CSF1R, 0.75),      # Third quartile (Q3)
         IQR = Q3 - Q1,                   # Calculate IQR
         Lower_Bound = Q1 - 1.5 * IQR,    # Lower bound
         Upper_Bound = Q3 + 1.5 * IQR) %>%
  filter(CSF1R >= Lower_Bound & CSF1R <= Upper_Bound) %>%  # Filter out outliers
  select(-Q1, -Q3, -IQR, -Lower_Bound, -Upper_Bound)       # Remove auxiliary columns

df_CSF1R <- as.data.frame(df_CSF1R)
df_CSF1R %>% 
  group_by(diagnosis) %>% 
  get_summary_stats(CSF1R, type = "mean_sd")
leveneTest(CSF1R ~ diagnosis, data = df_CSF1R)

bxp_CSF1R <- ggboxplot(
  df_CSF1R, x = "diagnosis", y = "CSF1R", 
  ylab = "CSF1R", xlab = "Diagnosis", add = "jitter")
bxp_CSF1R

stat.test_CSF1R <- df_CSF1R %>% 
  t_test(CSF1R ~ diagnosis, var.equal = TRUE) %>%
  add_significance()
stat.test_CSF1R

stat.test_CSF1R <- stat.test_CSF1R %>% add_xy_position(x = "diagnosis")
x_CSF1R  <- bxp_CSF1R + 
  stat_pvalue_manual(stat.test_CSF1R, tip.length = 0) +
  labs(subtitle = get_test_label(stat.test_CSF1R, detailed = TRUE))
x_CSF1R 
ggsave("ttest_CSF1R.png", plot = x_CSF1R, width = 12, height = 10, dpi = 300)

























