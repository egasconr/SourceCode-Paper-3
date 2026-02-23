# ==============================================================================
# Script: Venn Diagram Analysis
# ==============================================================================

#Library
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")

library(tidyverse)
library("ggVennDiagram")
library(famSKATRC)
library(VennDiagram)
library(gplots)
library(ggvenn)

# Open all files and filter

# --- Alzheimer's Disease (AD) ---
AD <- read.csv("AD_DEG.csv", header = TRUE, sep=",")
colnames(AD)[4] <- "SYMBOL"
AD <- AD %>% 
  select(-c(X, GENE.DESCRIPTION)) %>% 
  relocate("gene_type") %>% 
  relocate("SYMBOL")

# --- Amyotrophic Lateral Sclerosis (ALS) ---
ALS <- read.csv("ALS_DEG.csv", header = TRUE, sep=";")
colnames(ALS)[2] <- "SYMBOL"
ALS <- ALS %>% 
  select(-c(X)) %>% 
  relocate("gene_type") %>% 
  relocate("SYMBOL")


# --- Parkinson's Disease (PD) ---
PD <- read.csv("PD_DEG.csv", header = TRUE, sep=";")
PD <- PD %>% 
  select(-c(X, ID)) %>% 
  relocate("gene_type") %>% 
  relocate("SYMBOL")

# Filter and create a list selecting all Up and Down genes for each disease
AD_filter <- AD %>% 
  filter(AD$gene_type == "Up" | AD$gene_type == "Down")
AD_list <- AD_filter$SYMBOL

ALS_filter <- ALS %>% 
  filter(ALS$gene_type == "Up" | ALS$gene_type == "Down")
ALS_list <- ALS_filter$SYMBOL

PD_filter <- PD %>% 
  filter(PD$gene_type == "Up" | PD$gene_type == "Down")
PD_list <- PD_filter$SYMBOL

# Group everything into three groups
x <- list(
  AD = AD_list,
  ALS = ALS_list,
  PD = PD_list)

# Representation

venn.diagram(x,
             category.names = c("Alzheimer" , "Sclerosis Lateral Amyotrophic" , "Parkinson"),
             filename = 'venn.tiff',
             output = TRUE ,
             imagetype="tiff" ,
             height = 1000 , 
             width = 1000 , 
             resolution = 300,
             compression = "lzw",
             lwd = 1,
             col=c("#440154ff", '#21908dff', '#fde725ff'),
             fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
             cex = 0.5,
             fontfamily = "sans",
             cat.cex = 0.8,
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
             rotation = 1
)


# Extract common elements of the three groups (VennDiagram and gplots Packages):

ItemsList <- venn(x, show.plot = FALSE)
lenghts <- lengths(attributes(ItemsList)$intersections)
# lenghts

attributes <- attributes(ItemsList)$intersections
sink("Venn_info.txt")
print(lenghts)
print(attributes)
sink()
