
# Changes directory ----
#n = "/Users/Dao/Documents/MIT Grad/Classes/20.440/Project"
#setwd(dn)

# Loads in the Data----
library(Seurat)
library(DESeq2)
library(magrittr)
library(dplyr)
library(ggplot2)
library(escape)
library(SingleCellExperiment)
library(dittoSeq)
library(SeuratObject)

ts5 <- read.csv('Data/aaq0681_TableS5.csv') #Healthy
ts7 <- read.csv('Data/aaq0681_TableS7.csv') # Connective Tissue
ts9 <- read.csv('Data/Table_S9.csv') #Post amputation


# Processes the data ----
#Combine datasets
raw_counts <- rbind(ts5,ts9);
ts5_sub <- ts5[, !names(ts5) %in% 
                 c('cell_id','ident', 'orig.ident',
                   'tSNE_1', 'tSNE_2', 'nGene')]
ts7_sub <- ts7[, !names(ts7) %in% c('cell_id', 'celltype', 'timepoint')]
rownames(ts5_sub) <- ts5[,1]
rownames(ts7_sub) <- ts7[,1]
rownames(ts9_sub) <- ts9[,1]

ts5_SO <- CreateSeuratObject(counts = t(ts5_sub), min.cells = 3, min.features  = 100)
ts7_SO <- CreateSeuratObject(counts = t(ts7_sub), min.cells = 3, min.features  = 100)
ts9_SO <- CreateSeuratObject(counts = t(ts9_sub), min.cells = 3, min.features  = 100)

merged_CT <- merge(ts5_SO, y=ts7_SO)

#Creates meta data table
meta_data <- raw_counts[1:6]
colnames(meta_data) <- colnames(raw_counts)[1:6]

#Creates counts table (gene x library/cell)
count_table <- raw_counts[, !names(raw_counts) %in% 
                            c('cell_id','ident', 'orig.ident',
                              'tSNE_1', 'tSNE_2', 'nGene')]

#raw_counts <- raw_counts[oncs$Oncogene, ]
rownames(count_table) <- raw_counts[,1]
count_table <- t(count_table)

#count_table <- count_table[rownames(count_table) %in% oncs$Oncogene]

# Creates Seurat Object ----
# https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
SO <- CreateSeuratObject(counts = count_table, min.cells = 3, min.features  = 100)
SO <- NormalizeData(SO)
SO@meta.data['orig.ident'] <- meta_data[,3]
SO <- FindVariableFeatures(object = SO, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
# head(x = HVFInfo(object = SO)) #View the variance
SO <- ScaleData(SO, features = rownames(SO))

# Removes cell-cycle variation ----
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Gets cell cycle score and changes the current identity
SO <- CellCycleScoring(SO, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(SO, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# "Before" PCA
SO <- RunPCA(SO, features = c(s.genes, g2m.genes))
DimPlot(SO)

# Regressed PCA
SO <- ScaleData(SO, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SO))
SO <- RunPCA(SO, features = c(s.genes, g2m.genes))
DimPlot(SO)

# Make plot of all cells
SO <- RunUMAP(SO)

# Define Oncogene subset
oncs <- read.table("Data/axolotl_oncogenes.csv", header=T)
oncs <- oncs$Oncogene
SO_oncs <- subset(SO, features = oncs)
SO_oncs <- RunUMAP(SO_oncs)

# Make plot of uninjured cells
g1_treat <- WhichCells(SO_oncs, idents = c("CTa", "CTb"))
g1_untreat <- g1_treat <- WhichCells(SO_oncs, idents = c("18dpaA", "18dpaB", "25dpaA", "25dpaB", "38dpaB", "38dpaA"))
DimPlot(SO_oncs, label=T, group.by="identity", cells.highlight= list(g1_treat, g1_untreat), cols.highlight = c("darkblue", "darkred"), cols= "grey")

# Make plot of 18dpa cells
g1_treat <- WhichCells(SO_oncs, idents = c("18dpaA", "18dpaB"))
g1_untreat <- g1_treat <- WhichCells(SO_oncs, idents = c("CTa", "CTb", "25dpaA", "25dpaB", "38dpaB", "38dpaA"))
DimPlot(SO_oncs, label=T, group.by="identity", cells.highlight= list(g1_treat, g1_untreat), cols.highlight = c("darkblue", "darkred"), cols= "grey")

# Make plot of 25dpa cells
g1_treat <- WhichCells(SO_oncs, idents = c("25dpaA", "25dpaB"))
g1_untreat <- g1_treat <- WhichCells(SO_oncs, idents = c("CTa", "CTb", "18dpaA", "18dpaB", "38dpaB", "38dpaA"))
DimPlot(SO_oncs, label=T, group.by="identity", cells.highlight= list(g1_treat, g1_untreat), cols.highlight = c("darkblue", "darkred"), cols= "grey")

# Make plot of 38dpa cells
g1_treat <- WhichCells(SO_oncs, idents = c("38dpaA", "38dpaB"))
g1_untreat <- g1_treat <- WhichCells(SO_oncs, idents = c("CTa", "CTb", "25dpaA", "25dpaB", "18dpaB", "18dpaA"))
DimPlot(SO_oncs, label=T, group.by="identity", cells.highlight= list(g1_treat, g1_untreat), cols.highlight = c("darkblue", "darkred"), cols= "grey")

# Fig S2
# https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
merged_CT <- NormalizeData(merged_CT)
# head(x = HVFInfo(object = SO)) #View the variance
merged_CT <- ScaleData(merged_CT, features = rownames(SO))

# Removes cell-cycle variation ----
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Gets cell cycle score and changes the current identity
merged_CT <- CellCycleScoring(merged_CT, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(SO, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# "Before" PCA
merged_CT <- RunPCA(merged_CT, features = c(s.genes, g2m.genes))
DimPlot(merged_CT)

# Regressed PCA
merged_CT <- ScaleData(merged_CT, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SO))
merged_CT <- RunPCA(merged_CT, features = c(s.genes, g2m.genes))
DimPlot(merged_CT)

# Make plot of all cells
merged_CT <- RunUMAP(merged_CT)
FeaturePlot(merged_CT, features=c("FZD2", "KLF4", "FZD7"))

# Rename identities
SO$orig.ident[SO$orig.ident %in% c("CTb", "CTa")] = "Ctrl"
SO$orig.ident[SO$orig.ident %in% c("38dpaA", "38dpaB")] = "38dpa"
SO$orig.ident[SO$orig.ident %in% c("25dpaA", "25dpaB")] = "25dpa"
SO$orig.ident[SO$orig.ident %in% c("18dpaA", "18dpaB")] = "18dpa"

# Do single single gsea
GS.hallmark <- getGeneSets(library = "H")
ES.seurat <- enrichIt(obj = SO, gene.sets = GS.hallmark, groups = 1000, cores = 2)

SO <- Seurat::AddMetaData(SO, ES.seurat)
dittoHeatmap(SO, genes = NULL, metas = names(ES.seurat))
