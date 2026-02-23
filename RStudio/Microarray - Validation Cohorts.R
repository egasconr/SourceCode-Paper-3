# ==============================================================================
# Script: Complete Analysis of Validation Cohorts (AD, ALS, PD)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Package Loading
# ------------------------------------------------------------------------------

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
library(gridExtra)
library(pROC)
library(limma)
library(tibble)
library(ggpubr)
library(rstatix)
library(car)
library(sva) # Required for ComBat

# ==============================================================================
# 2. ALZHEIMER'S DISEASE (AD) COHORT - GSE63061
# ==============================================================================

# ------------------------------------------------------------------------------
# 2.1. Metadata Processing
# ------------------------------------------------------------------------------

# Load GEO dataset
GSE63061 <- getGEO('GSE63061', GSEMatrix = TRUE)
length(GSE63061)
class(GSE63061[[1]])

# Prepare and analyze metadata labels
varLabels(GSE63061[[1]]) 
info_AD <- pData(GSE63061[[1]]) 

# Select relevant columns (Status:ch1)
# write.csv(info_AD, file = "info_AD_metadata.csv")
info_AD <- info_AD[,c(1,2,39)]
colnames(info_AD)

info_AD <- info_AD %>%
  rename(diagnosis = "status:ch1")

# ------------------------------------------------------------------------------
# 2.2. Expression Data Preprocessing
# ------------------------------------------------------------------------------

# Save expression data to object
gse_AD <- GSE63061[[1]]

# Analyze expression data distribution
summary(exprs(gse_AD))[,1:200]
boxplot(exprs(gse_AD)[,1:200], outline=FALSE, ylab = "Expression level", xlab = "Samples", col = "red")
abline(h=median(exprs(gse_AD)),col="blue")

# Create dataframe
gse.expr_AD <- exprs(gse_AD)
normalized.expr.dataframe_AD <- as.data.frame(gse.expr_AD)

# Save original dataframe
# write.csv(normalized.expr.dataframe_AD, file = "normalized.exprs.dataframe_AD.csv")

# Annotate probes with SYMBOL
feature.data_AD <- GSE63061$GSE63061_series_matrix.txt.gz@featureData@data
feature.data.selection_AD <- feature.data_AD[,c(1,13)]
colnames(feature.data.selection_AD)[2] <- "SYMBOL"

# Merge dataframes
normalized.expr.dataframe.ann_AD <- normalized.expr.dataframe_AD %>%
  rownames_to_column(var = "ID") %>%
  inner_join(.,feature.data.selection_AD, by = "ID") 

# Check structure
colnames(normalized.expr.dataframe.ann_AD)
# write.csv(normalized.expr.dataframe.ann_AD, "normalized.expr.annotated_AD.csv")

# ------------------------------------------------------------------------------
# 2.3. Filtering and Cleaning
# ------------------------------------------------------------------------------

# Modify dataframe to have SYMBOL as ID and remove probes
data_AD <- normalized.expr.dataframe.ann_AD %>% 
  select(-c(ID)) %>% 
  relocate("SYMBOL") %>% 
  na.omit() 

# Remove duplicate genes
head(data_AD)[1:10]
has_rownames(data_AD)
data.mod_AD <- remove_rownames(data_AD) 
has_rownames(data.mod_AD)

# Dplyr requires unique identifiers
duplicado_AD <- duplicated(data.mod_AD$SYMBOL)
sum(duplicado_AD)

data.mod.uniq_AD <- data.mod_AD %>% 
  distinct(SYMBOL, .keep_all = TRUE)
sum(duplicated(data.mod.uniq_AD$SYMBOL))

data.mod.uniq_AD <- column_to_rownames(data.mod.uniq_AD, var = 'SYMBOL')

# Save unique data
# write.csv(data.mod.uniq_AD, "data.mod.uniq_AD.csv")

# Transpose and filter for specific conditions (AD vs CTL, removing MCI)
data.mod.t_AD <- t(data.mod.uniq_AD)

data.mod.t.filter_AD <- as.data.frame(cbind(data.mod.t_AD, info_AD[3])) # Add diagnosis info
data.mod.t.filter_AD <- data.mod.t.filter_AD %>% relocate("diagnosis")

# Filter samples
data.mod.t.filter_AD <- data.mod.t.filter_AD %>% 
  filter(data.mod.t.filter_AD$diagnosis == 'CTL' | data.mod.t.filter_AD$diagnosis == 'AD')

ncol(data.mod.t.filter_AD)
nrow(data.mod.t.filter_AD)

# write.csv(data.mod.t.filter_AD, "data.mod.t.filter_AD.csv")

# ------------------------------------------------------------------------------
# 2.4. Differential Expression Analysis (DEG) - AD
# ------------------------------------------------------------------------------

# Check patient counts
unique(gse_AD$`status:ch1`)
AD <- gse_AD$`status:ch1` == "AD"
CONTROL_AD <- gse_AD$`status:ch1` == "CTL"
sum(AD) 
sum(CONTROL_AD) 

# Design Matrix
experimental.design_AD <- model.matrix( ~ 0 + gse_AD[['status:ch1']])
colnames(experimental.design_AD) <- levels(as.factor(gse_AD[['status:ch1']]))
colnames(experimental.design_AD) <- c("AD", "borderline_MCI", "CONTROL", "CTL_to_AD", "MCI", "MCI_to_CTL", "OTHER") 
experimental.design_AD

# Contrast Matrix
contrast_matrix_AD <- makeContrasts(AD - CONTROL, levels=experimental.design_AD)
contrast_matrix_AD

# Compensate outliers with empirical weights
aw_AD <- arrayWeights(exprs(gse_AD), experimental.design_AD)

# Linear Fit
linear.fit_AD <- lmFit(gse_AD,experimental.design_AD, weights = aw_AD)
contrast.linear.fit_AD <- contrasts.fit(linear.fit_AD,contrasts=contrast_matrix_AD)
contrast.results_AD <- eBayes(contrast.linear.fit_AD)
summary(decideTests(contrast.results_AD,lfc=0.2))
topTable(contrast.results_AD)

# Extract Condition Data
Condition_Final_AD <- topTable(contrast.results_AD, number = 32049, coef = 1, sort.by = "logFC" )
head(Condition_Final_AD)

# Annotation and Cleaning
Annotated_Condition_Final_AD <- Condition_Final_AD[,c(6,31,32,33,34,35,36)]

# Clean empty gene symbols
Annotated_Condition_Final_AD$ILMN_Gene[Annotated_Condition_Final_AD$ILMN_Gene == " "] <- NA
colSums(is.na(Annotated_Condition_Final_AD))

# Remove duplicates in results
head(Annotated_Condition_Final_AD)
has_rownames(Annotated_Condition_Final_AD)
Annotated_Condition_Final_AD.mod <- remove_rownames(Annotated_Condition_Final_AD)
has_rownames(Annotated_Condition_Final_AD.mod)

duplicado_AD <- duplicated(Annotated_Condition_Final_AD.mod$ILMN_Gene)
sum(duplicado_AD)

Annotated_Condition_Final_AD.mod <- Annotated_Condition_Final_AD.mod %>% 
  distinct(ILMN_Gene, .keep_all = TRUE)
sum(duplicated(Annotated_Condition_Final_AD.mod$ILMN_Gene))

write.csv(Annotated_Condition_Final_AD.mod, "Annotated_Condition_Final_AD.csv")

# ------------------------------------------------------------------------------
# 2.5. Visualization - AD
# ------------------------------------------------------------------------------

# Volcano Plot Preview
volplot <- ggplot(Annotated_Condition_Final_AD.mod, aes(x = logFC, y= -log10(adj.P.Val))) + geom_point()
volplot +  geom_hline(yintercept = -log10(0.05),linetype = "dashed") + 
  geom_vline(xintercept = c(log2(1.07), log2(0.93)),linetype = "dashed") + 
  xlim(-1, 1)

# Categorize genes (Up/Down/NS) with Cutoff 0.1
Annotated_Condition_Final_DEG_AD <- Annotated_Condition_Final_AD.mod %>%
  mutate(gene_type = case_when(Annotated_Condition_Final_AD.mod$logFC > 0.1 & Annotated_Condition_Final_AD.mod$adj.P.Val < 0.05 ~ "Up",
                               Annotated_Condition_Final_AD.mod$logFC < -0.1 & Annotated_Condition_Final_AD.mod$adj.P.Val < 0.05 ~ "Down",
                               TRUE ~ "NS"))   

# Counts
Annotated_Condition_Final_DEG_AD %>%
  count(gene_type)

# Save Results
write.csv(Annotated_Condition_Final_DEG_AD, "Validation_AD_DEG.csv")
write.table(Annotated_Condition_Final_DEG_AD, "Validation_AD_DEG.csv", sep = ";", dec = ",", row.names = TRUE)

# Final Volcano Plot
cols <- c("Up" = "red", "Down" = "blue", "NS" = "grey") 
sizes <- c("Up" = 3, "Down" =3, "NS" = 2) 
alphas <- c("Up" = 1, "Down" = 1, "NS" = 0.5)

Annotated_Condition_Final_DEG_AD %>%
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
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) + 
  scale_alpha_manual(values = alphas) + xlim(-1, 1) + ylim(0,16) + 
  labs(title = " ", 
       x = "log2(fold change)",
       y = "−log10(adjusted P.value)",
       font.size = 16,
       colour = "Expression \nchange") +
  theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

# ==============================================================================
# 3. ALS COHORT - GSE112680
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1. Metadata Processing
# ------------------------------------------------------------------------------

# Load GEO dataset
GSE112680 <- getGEO('GSE112680', GSEMatrix = TRUE)
length(GSE112680)
class(GSE112680[[1]])

# Prepare metadata
varLabels(GSE112680[[1]]) 
info_ALS <- pData(GSE112680[[1]]) 

# Select columns (Diagnosis)
# write.csv(info_ALS, file = "info_ALS_metadata.csv")
info_ALS <- info_ALS[,c(1,2,42)]
colnames(info_ALS)

info_ALS <- info_ALS %>%
  rename(diagnosis = "diagnosis:ch1")

# ------------------------------------------------------------------------------
# 3.2. Expression Data Preprocessing
# ------------------------------------------------------------------------------

gse_ALS <- GSE112680[[1]]

# Analyze distribution
summary(exprs(gse_ALS))[,1:200] 
boxplot(exprs(gse_ALS)[,1:200], outline=FALSE, ylab = "Expression level", xlab = "Samples", col = "red")
abline(h=median(exprs(gse_ALS)),col="blue")

gse.expr_ALS <- exprs(gse_ALS)
normalized.expr.dataframe_ALS <- as.data.frame(gse.expr_ALS)

# write.csv(normalized.expr.dataframe_ALS, file = "normalized.exprs.dataframe_ALS.csv")

# Annotate probes with SYMBOL
feature.data_ALS <- GSE112680$GSE112680_series_matrix.txt.gz@featureData@data
feature.data.selection_ALS <- feature.data_ALS[,c(1,13)]
colnames(feature.data.selection_ALS)[2] <- "SYMBOL"

# Merge
normalized.expr.dataframe.ann_ALS <- normalized.expr.dataframe_ALS %>%
  rownames_to_column(var = "ID") %>%
  inner_join(.,feature.data.selection_ALS, by = "ID") 

# Check
colnames(normalized.expr.dataframe.ann_ALS)
# write.csv(normalized.expr.dataframe.ann_ALS, "normalized.expr.annotated_ALS.csv")

# ------------------------------------------------------------------------------
# 3.3. Filtering and Cleaning
# ------------------------------------------------------------------------------

data_ALS <- normalized.expr.dataframe.ann_ALS %>% 
  select(-c(ID)) %>% 
  relocate("SYMBOL") %>% 
  na.omit() 

# Remove duplicates
head(data_ALS)[1:10]
has_rownames(data_ALS)
data.mod_ALS <- remove_rownames(data_ALS) 
has_rownames(data.mod_ALS)

duplicado_ALS <- duplicated(data.mod_ALS$SYMBOL)
sum(duplicado_ALS)

data.mod.uniq_ALS <- data.mod_ALS %>% 
  distinct(SYMBOL, .keep_all = TRUE)
sum(duplicated(data.mod.uniq_ALS$SYMBOL))

data.mod.uniq_ALS <- column_to_rownames(data.mod.uniq_ALS, var = 'SYMBOL')

# write.csv(data.mod.uniq_ALS, "data.mod.uniq_ALS.csv")

# Transpose and Filter samples (Remove MIM, keep CON and ALS)
data.mod.t_ALS <- t(data.mod.uniq_ALS)

data.mod.t.filter_ALS <- as.data.frame(cbind(data.mod.t_ALS, info_ALS[3])) 
data.mod.t.filter_ALS <- data.mod.t.filter_ALS %>% relocate("diagnosis")

data.mod.t.filter_ALS <- data.mod.t.filter_ALS %>% 
  filter(data.mod.t.filter_ALS$diagnosis == 'CON' | data.mod.t.filter_ALS$diagnosis == 'ALS')

ncol(data.mod.t.filter_ALS)
nrow(data.mod.t.filter_ALS)

# write.csv(data.mod.t.filter_ALS, "data.mod.t.filter_ALS.csv")

# ------------------------------------------------------------------------------
# 3.4. Differential Expression Analysis (DEG) - ALS
# ------------------------------------------------------------------------------

# Check counts
unique(gse_ALS$`diagnosis:ch1`)
ALS <- gse_ALS$`diagnosis:ch1` == "ALS"
CONTROL_ALS <- gse_ALS$`diagnosis:ch1` == "CON"
sum(ALS) 
sum(CONTROL_ALS) 

# Design Matrix
experimental.design_ALS <- model.matrix( ~ 0 + gse_ALS[['diagnosis:ch1']])
colnames(experimental.design_ALS) <- levels(as.factor(gse_ALS[['diagnosis:ch1']]))
experimental.design_ALS

# Contrast Matrix
contrast_matrix_ALS <- makeContrasts(ALS - CON, levels=experimental.design_ALS)
contrast_matrix_ALS

# Compensate outliers
aw_ALS <- arrayWeights(exprs(gse_ALS), experimental.design_ALS)

# Linear Fit
linear.fit_ALS <- lmFit(gse_ALS,experimental.design_ALS, weights = aw_ALS)
contrast.linear.fit_ALS <- contrasts.fit(linear.fit_ALS,contrasts=contrast_matrix_ALS)
contrast.results_ALS <- eBayes(contrast.linear.fit_ALS)
summary(decideTests(contrast.results_ALS,lfc=0.2))
topTable(contrast.results_ALS)

# Extract Condition Data
Condition_Final_ALS <- topTable(contrast.results_ALS, number = 29830, coef = 1, sort.by = "logFC" )
head(Condition_Final_ALS)

# Annotation and Cleaning
Annotated_Condition_Final_ALS <- Condition_Final_ALS[,c(6,31,32,33,34,35,36)]

Annotated_Condition_Final_ALS$ILMN_Gene[Annotated_Condition_Final_ALS$ILMN_Gene == " "] <- NA
colSums(is.na(Annotated_Condition_Final_ALS))

# Remove duplicates
head(Annotated_Condition_Final_ALS)
has_rownames(Annotated_Condition_Final_ALS)
Annotated_Condition_Final_ALS.mod <- remove_rownames(Annotated_Condition_Final_ALS)
has_rownames(Annotated_Condition_Final_ALS.mod)

duplicado_ALS <- duplicated(Annotated_Condition_Final_ALS.mod$ILMN_Gene)
sum(duplicado_ALS)

Annotated_Condition_Final_ALS.mod <- Annotated_Condition_Final_ALS.mod %>% 
  distinct(ILMN_Gene, .keep_all = TRUE)
sum(duplicated(Annotated_Condition_Final_ALS.mod$ILMN_Gene))

write.csv(Annotated_Condition_Final_ALS.mod, "Annotated_Condition_Final_ALS.csv")
write.table(Annotated_Condition_Final_ALS.mod, "Annotated_Condition_Final_ALS.csv", sep = ";", dec = ",", row.names = TRUE)

# ------------------------------------------------------------------------------
# 3.5. Visualization - ALS
# ------------------------------------------------------------------------------

# Volcano Plot Preview
volplot <- ggplot(Annotated_Condition_Final_ALS.mod, aes(x = logFC, y= -log10(adj.P.Val))) + geom_point()
volplot +  geom_hline(yintercept = -log10(0.05),linetype = "dashed") + 
  geom_vline(xintercept = c(log2(1.148), log2(0.87)),linetype = "dashed") + 
  xlim(-1, 1)

# Categorize genes (Up/Down/NS) with Cutoff 0.2
Annotated_Condition_Final_DEG_ALS <- Annotated_Condition_Final_ALS.mod %>%
  mutate(gene_type = case_when(Annotated_Condition_Final_ALS.mod$logFC > 0.2 & Annotated_Condition_Final_ALS.mod$adj.P.Val < 0.05 ~ "Up",
                               Annotated_Condition_Final_ALS.mod$logFC < -0.2 & Annotated_Condition_Final_ALS.mod$adj.P.Val < 0.05 ~ "Down",
                               TRUE ~ "NS"))   

# Counts
Annotated_Condition_Final_DEG_ALS %>%
  count(gene_type)

# Save Results
write.csv(Annotated_Condition_Final_DEG_ALS, "Validation_ALS_DEG.csv")
write.table(Annotated_Condition_Final_DEG_ALS, "Validation_ALS_DEG.csv", sep = ";", dec = ",", row.names = TRUE)

# Final Volcano Plot
cols <- c("Up" = "red", "Down" = "blue", "NS" = "grey") 
sizes <- c("Up" = 3, "Down" =3, "NS" = 2) 
alphas <- c("Up" = 1, "Down" = 1, "NS" = 0.5)

Annotated_Condition_Final_DEG_ALS %>%
  ggplot(aes(x = logFC, 
             y= -log10(adj.P.Val),
             fill = gene_type,    
             size = gene_type,
             alpha = gene_type)) + 
  geom_point(shape = 21, 
             colour = "black") + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.87), log2(1.148)),
             linetype = "dashed") +
  scale_fill_manual(values = cols) + 
  scale_size_manual(values = sizes) + 
  scale_alpha_manual(values = alphas) + xlim(-1, 1) + ylim(0,25) + 
  labs(title = " ", 
       x = "log2(fold change)",
       y = "−log10(adjusted P.value)",
       font.size = 16,
       colour = "Expression \nchange") +
  theme_bw() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

# ==============================================================================
# 4. PARKINSON'S DISEASE (PD) COHORTS - GSE57475, GSE72267, GSE18838
# ==============================================================================

# ------------------------------------------------------------------------------
# 4.1. Metadata Processing
# ------------------------------------------------------------------------------

# --- PD1: GSE57475 ---
GSE57475 <- getGEO('GSE57475', GSEMatrix = TRUE)
info_PD1 <- pData(GSE57475[[1]])
# write.csv(info_PD1, file = "info_PD1_metadata.csv")
info_PD1 <- info_PD1[,c(1,2,34)]
colnames(info_PD1)
info_PD1 <- info_PD1 %>% rename(diagnosis = "disease state:ch1")

# --- PD2: GSE72267 ---
GSE72267 <- getGEO('GSE72267', GSEMatrix = TRUE)
info_PD2 <- pData(GSE72267[[1]])
# write.csv(info_PD2, file = "info_PD2_metadata.csv")
info_PD2 <- info_PD2[,c(1,2,32)]
colnames(info_PD2)
info_PD2 <- info_PD2 %>% rename(diagnosis = "diagnosis:ch1")

# --- PD3: GSE18838 ---
GSE18838 <- getGEO('GSE18838', GSEMatrix = TRUE)
info_PD3 <- pData(GSE18838[[1]])
# write.csv(info_PD3, file = "info_PD3_metadata.csv")
info_PD3 <- info_PD3[,c(1,2,35)]
colnames(info_PD3)
info_PD3 <- info_PD3 %>% rename(diagnosis = "disease state:ch1")

# Combine metadata (if needed for reference, but we process separately first)
info_PD_total <- rbind(info_PD1, info_PD2, info_PD3)

# ------------------------------------------------------------------------------
# 4.2. Expression Data Preprocessing per Dataset
# ------------------------------------------------------------------------------

# --- PD1: GSE57475 (From .txt matrix) ---
gse_PD1 <- GSE57475[[1]]

# Normalization (Quantile normalization on Log2 transformed data)
gse.expr_PD1 <- log2(exprs(gse_PD1)) 
normalized.gse_PD1 <- normalizeBetweenArrays(gse.expr_PD1, method="quantile")
normalized.expr.dataframe_PD1 <- as.data.frame(normalized.gse_PD1)
# write.csv(normalized.expr.dataframe_PD1, file = "normalized.exprs.dataframe_PD1.csv") 

# Annotation
feature.data_PD1 <- GSE57475$GSE57475_series_matrix.txt.gz@featureData@data
feature.data.selection_PD1 <- feature.data_PD1[,c(1,14)]
colnames(feature.data.selection_PD1)[2] <- "SYMBOL"

normalized.expr.dataframe.ann_PD1 <- normalized.expr.dataframe_PD1 %>%
  rownames_to_column(var = "ID") %>%
  inner_join(.,feature.data.selection_PD1, by = "ID") 

# Cleaning PD1
data_PD1 <- normalized.expr.dataframe.ann_PD1 %>% 
  select(-c(ID)) %>% 
  relocate("SYMBOL") %>% 
  na.omit() 

data.mod_PD1 <- remove_rownames(data_PD1) 
data.mod.uniq_PD1 <- data.mod_PD1 %>% 
  distinct(SYMBOL, .keep_all = TRUE) %>%
  column_to_rownames(var = 'SYMBOL')
# write.csv(data.mod.uniq_PD1, "data.mod.uniq_PD1.csv")

# Filter PD1 Samples
data.mod.t_PD1 <- t(data.mod.uniq_PD1)
data.mod.t.filter_PD1 <- as.data.frame(cbind(data.mod.t_PD1, info_PD1[3])) 
data.mod.t.filter_PD1 <- data.mod.t.filter_PD1 %>% relocate("diagnosis")
# write.csv(data.mod.t.filter_PD1, "data.mod.t.filter_PD1.csv")


# --- PD2: GSE72267 (From .CEL files) ---
# Note: Ensure RAW data is downloaded
GSE72267_obj <- GSE72267[[1]]
getGEOSuppFiles("GSE72267")
untar("GSE72267/GSE72267_RAW.tar", exdir = 'data/')
gse_PD2 <- ReadAffy(celfile.path = "data/", phenoData=phenoData(GSE72267_obj))

# Normalization (RMA)
normalized.gse_PD2 <- affy::rma(gse_PD2)
normalized.expr.dataframe_PD2 <- as.data.frame(exprs(normalized.gse_PD2))
# write.csv(normalized.expr.dataframe_PD2, file = "normalized.exprs.dataframe_PD2.csv") 

# Annotation
feature.data_PD2 <- GSE72267$GSE72267_series_matrix.txt.gz@featureData@data
feature.data.selection_PD2 <- feature.data_PD2[,c(1,11)]
colnames(feature.data.selection_PD2)[2] <- "SYMBOL"
feature.data.selection_PD2 <- feature.data.selection_PD2 %>%
  mutate(SYMBOL = gsub("///.*", "", SYMBOL))

normalized.expr.dataframe.ann_PD2 <- normalized.expr.dataframe_PD2 %>%
  rownames_to_column(var = "ID") %>%
  inner_join(.,feature.data.selection_PD2, by = "ID") 

# Cleaning PD2
data_PD2 <- normalized.expr.dataframe.ann_PD2 %>% 
  select(-c(ID)) %>% 
  relocate("SYMBOL") %>% 
  na.omit() 

data.mod_PD2 <- remove_rownames(data_PD2) 
data.mod.uniq_PD2 <- data.mod_PD2 %>% 
  distinct(SYMBOL, .keep_all = TRUE) %>%
  column_to_rownames(var = 'SYMBOL')
# write.csv(data.mod.uniq_PD2, "data.mod.uniq_PD2.csv")

# Filter PD2 Samples
data.mod.t_PD2 <- t(data.mod.uniq_PD2)
data.mod.t.filter_PD2 <- as.data.frame(cbind(data.mod.t_PD2, info_PD2[3])) 
data.mod.t.filter_PD2 <- data.mod.t.filter_PD2 %>% relocate("diagnosis")
# write.csv(data.mod.t.filter_PD2, "data.mod.t.filter_PD2.csv")


# --- PD3: GSE18838 (From .CEL files) ---
# Note: Ensure RAW data is downloaded
GSE18838_obj <- GSE18838[[1]]
getGEOSuppFiles("GSE18838")
untar("GSE18838/GSE18838_RAW.tar", exdir = 'dataPD3/')

# Load .CEL files (using read.celfiles for specific platform compatibility)
cel_files_PD3 <- list.celfiles("dataPD3/", full.names = TRUE)
gse_PD3 <- read.celfiles(cel_files_PD3)

# Associate Metadata
pheno_data <- pData(phenoData(GSE18838_obj))
rownames(pheno_data) <- sampleNames(gse_PD3)
pData(gse_PD3) <- pheno_data

# Normalization (RMA via Oligo)
normalized.gse_PD3 <- oligo::rma(gse_PD3)
normalized.expr.dataframe_PD3 <- as.data.frame(exprs(normalized.gse_PD3))
# write.csv(normalized.expr.dataframe_PD3, file = "normalized.exprs.dataframe_PD3.csv") 

# Annotation
feature.data_PD3 <- GSE18838$GSE18838_series_matrix.txt.gz@featureData@data
feature.data.selection_PD3 <- feature.data_PD3[,c(1,10)]
colnames(feature.data.selection_PD3)[2] <- "SYMBOL"

# Clean Symbol splitting " // "
feature.data.selection_PD3 <- feature.data.selection_PD3 %>%
  mutate(SYMBOL = sapply(str_split(SYMBOL, " // "), function(x) x[2])) %>% 
  mutate(ID = as.character(ID))

normalized.expr.dataframe.ann_PD3 <- normalized.expr.dataframe_PD3 %>%
  rownames_to_column(var = "ID") %>%
  inner_join(.,feature.data.selection_PD3, by = "ID") 

# Cleaning PD3
data_PD3 <- normalized.expr.dataframe.ann_PD3 %>% 
  select(-c(ID)) %>% 
  relocate("SYMBOL") %>% 
  na.omit() 

data.mod_PD3 <- remove_rownames(data_PD3) 
data.mod.uniq_PD3 <- data.mod_PD3 %>% 
  distinct(SYMBOL, .keep_all = TRUE) %>%
  column_to_rownames(var = 'SYMBOL')
# write.csv(data.mod.uniq_PD3, "data.mod.uniq_PD3.csv")

# Filter PD3 Samples
data.mod.t_PD3 <- t(data.mod.uniq_PD3)
data.mod.t.filter_PD3 <- as.data.frame(cbind(data.mod.t_PD3, info_PD3[3])) 
data.mod.t.filter_PD3 <- data.mod.t.filter_PD3 %>% relocate("diagnosis")
# write.csv(data.mod.t.filter_PD3, "data.mod.t.filter_PD3.csv")


# ------------------------------------------------------------------------------
# 4.3. Merging Datasets and Batch Correction
# ------------------------------------------------------------------------------

# Extract Common Genes
genes_PD1 <- colnames(data.mod.t.filter_PD1)[-which(colnames(data.mod.t.filter_PD1) == "diagnosis")]
genes_PD2 <- colnames(data.mod.t.filter_PD2)[-which(colnames(data.mod.t.filter_PD2) == "diagnosis")]
genes_PD3 <- colnames(data.mod.t.filter_PD3)[-which(colnames(data.mod.t.filter_PD3) == "diagnosis")]

common_genes <- Reduce(intersect, list(genes_PD1, genes_PD2, genes_PD3))

# Filter datasets for common genes
PD1_filtered <- data.mod.t.filter_PD1 %>% select(diagnosis, all_of(common_genes))
PD2_filtered <- data.mod.t.filter_PD2 %>% select(diagnosis, all_of(common_genes))
PD3_filtered <- data.mod.t.filter_PD3 %>% select(diagnosis, all_of(common_genes))

# Add Batch Information
PD1_filtered <- PD1_filtered %>% mutate(batch = "PD1")
PD2_filtered <- PD2_filtered %>% mutate(batch = "PD2")
PD3_filtered <- PD3_filtered %>% mutate(batch = "PD3")

# Standardize Diagnosis Labels & Ensure Numeric Data
# Note: gsub parentheses must be escaped
PD1_filtered <- PD1_filtered %>%
  mutate(diagnosis = gsub("Dopamine Transporter Imaging \\(DAT\\)-confirmed PD", "PD", diagnosis)) %>%
  mutate(diagnosis = gsub("healthy control", "CONTROL", diagnosis)) %>%
  as.data.frame()

PD2_filtered <- PD2_filtered %>%
  mutate(across(-c(diagnosis, batch), as.numeric)) %>% 
  mutate(diagnosis = gsub("Parkinson's disease", "PD", diagnosis)) %>%
  mutate(diagnosis = gsub("Healthy", "CONTROL", diagnosis)) %>%
  as.data.frame()

PD3_filtered <- PD3_filtered %>%
  mutate(across(-c(diagnosis, batch), as.numeric)) %>% 
  mutate(diagnosis = gsub("Parkinson's Disease", "PD", diagnosis)) %>%
  mutate(diagnosis = gsub("none", "CONTROL", diagnosis)) %>%
  as.data.frame()

# Combine Datasets
data_PD_total <- bind_rows(PD1_filtered, PD2_filtered, PD3_filtered)

# Separate diagnosis and batch vectors
diagnosis <- data_PD_total$diagnosis
batch <- data_PD_total$batch

# Expression Matrix for ComBat (Genes x Samples)
# Transpose because ComBat expects genomic features as rows
expr_matrix <- data_PD_total %>%
  select(-diagnosis, -batch) %>%
  as.matrix() %>%
  t() 

# Batch Correction using ComBat
mod <- model.matrix(~ as.factor(diagnosis))

corrected_expr_matrix <- ComBat(
  dat = expr_matrix,      
  batch = batch,          
  mod = mod,              
  par.prior = TRUE,       
  prior.plots = FALSE     
)

# Convert back to DataFrame (Samples x Genes)
corrected_expr_dataframe <- as.data.frame(t(corrected_expr_matrix))
corrected_expr_dataframe$diagnosis <- diagnosis

data.mod.t.filter_PD_total <- corrected_expr_dataframe %>% 
  relocate(diagnosis)

# Save combined corrected data
write.csv(data.mod.t.filter_PD_total, "PD_combined_corrected_data.csv", row.names = FALSE)

# PCA Visualization (Before vs After)
pca_before <- prcomp(t(expr_matrix), scale. = TRUE)
pca_after <- prcomp(t(corrected_expr_matrix), scale. = TRUE)

# Plot
par(mfrow=c(1,2))
plot(pca_before$x[, 1:2], col = as.factor(batch), main = "Parkinson's PCA before ComBat (Batch Color)")
plot(pca_after$x[, 1:2], col = as.factor(diagnosis), main = "Parkinson's PCA after ComBat (Diagnosis Color)")
par(mfrow=c(1,1))


# ==============================================================================
# 5. TARGETED EXPRESSION ANALYSIS (T-TESTS) - AD, ALS, PD
# ==============================================================================

# ------------------------------------------------------------------------------
# 5.1. AD Targeted Analysis
# ------------------------------------------------------------------------------

data.ttest_AD <- data.mod.t.filter_AD %>% 
  select(diagnosis, any_of(c("LILRB2", "TLN1", "ITGAM", "C5AR1", "RPL15", "FOS", "EIF4A2", "NKTR", "TPT1", "RBM25", "RPL22", "TAX1BP1", "RPS6KB1", "NCF4", "MYO1F", "VAV1", "HCLS1")))
write.csv(data.ttest_AD, "data.ttest_AD_validation.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 5.2. ALS Targeted Analysis
# ------------------------------------------------------------------------------

data.ttest_ALS <- data.mod.t.filter_ALS %>% 
  filter(diagnosis == "CON" | diagnosis == "ALS") %>% 
  select(diagnosis, any_of(c("LILRB2", "TLN1", "ITGAM", "C5AR1", "RPL15", "FOS", "EIF4A2", "NKTR", "TPT1", "RBM25", "RPL22", "TAX1BP1", "RPS6KB1", "IL18", "PTPRC", "PTEN", "CTSS", "TNFSF13B", "TLR4", "CSF1R")))
write.csv(data.ttest_ALS, "data.ttest_ALS_validation.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 5.3. PD Targeted Analysis (Example: MYO1F)
# ------------------------------------------------------------------------------

data.ttest_PD <- data.mod.t.filter_PD_total %>% 
  filter(diagnosis == "CONTROL" | diagnosis == "PD") %>% 
  select(diagnosis, any_of(c("LILRB2", "TLN1", "ITGAM", "C5AR1", "IL18", "PTPRC", "PTEN", "CTSS", "TLR4", "CSF1R", "NCF4", "MYO1F", "VAV1", "HCLS1")))
write.csv(data.ttest_PD, "data.ttest_PD_validation.csv", row.names = FALSE)

# Statistical Analysis for MYO1F
data.ttest_PD %>% 
  group_by(diagnosis) %>% 
  get_summary_stats(MYO1F, type = "mean_sd")

# Homogeneity of variance
leveneTest(MYO1F ~ diagnosis, data = data.ttest_PD) 

# Outlier Detection
data.ttest_PD %>%
  group_by(diagnosis) %>%
  identify_outliers(MYO1F) 

# Remove Outliers
df_MYO1F_PD <- data.ttest_PD %>%
  group_by(diagnosis) %>% 
  mutate(Q1 = quantile(MYO1F, 0.25),      
         Q3 = quantile(MYO1F, 0.75),      
         IQR = Q3 - Q1,                   
         Lower_Bound = Q1 - 1.5 * IQR,    
         Upper_Bound = Q3 + 1.5 * IQR) %>%
  filter(MYO1F >= Lower_Bound & MYO1F <= Upper_Bound) %>%  
  select(-Q1, -Q3, -IQR, -Lower_Bound, -Upper_Bound)       

df_MYO1F_PD <- as.data.frame(df_MYO1F_PD)

# Re-check stats
df_MYO1F_PD %>% 
  group_by(diagnosis) %>% 
  get_summary_stats(MYO1F, type = "mean_sd")

# T-test
stat.test_MYO1F_PD <- df_MYO1F_PD %>% 
  t_test(MYO1F ~ diagnosis, var.equal = TRUE) %>%
  add_significance()
stat.test_MYO1F_PD

# Boxplot
bxp_MYO1F <- ggboxplot(
  df_MYO1F_PD, x = "diagnosis", y = "MYO1F", 
  ylab = "MYO1F", xlab = "Diagnosis", add = "jitter")

# Final Plot with P-value
stat.test_MYO1F_PD <- stat.test_MYO1F_PD %>% add_xy_position(x = "diagnosis")
x_MYO1F  <- bxp_MYO1F + 
  stat_pvalue_manual(stat.test_MYO1F_PD, tip.length = 0) +
  labs(subtitle = get_test_label(stat.test_MYO1F_PD, detailed = TRUE))
x_MYO1F 

ggsave("ttest_MYO1F_PD.png", plot = x_MYO1F, width = 12, height = 10, dpi = 300)









