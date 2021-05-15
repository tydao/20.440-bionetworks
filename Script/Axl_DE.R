
# Changes directory ----
dn = "/Users/Dao/Documents/MIT Grad/Classes/20.440/Project"
setwd(dn)

# Loads in the Data----
library(Seurat)
library(SeuratDisk)
library(DESeq2)
library(magrittr)
library(dplyr)
library(ggplot2)

ts5 <- read.csv('Data/aaq0681_TableS5.csv') #Healthy
ts9 <- read.csv('Data/aa10681_TableS9.csv') #Post amputation


# Processes the data ----
#Combine datasets
raw_counts <- rbind(ts5,ts9);

#Creates meta data table
meta_data <- raw_counts[1:6]
colnames(meta_data) <- colnames(raw_counts)[1:6]

#Creates counts table (gene x library/cell)
count_table <- raw_counts[, !names(raw_counts) %in% 
                            c('cell_id','ident', 'orig.ident',
                              'tSNE_1', 'tSNE_2', 'nGene')]
rownames(count_table) <- raw_counts[,1]
count_table <- t(count_table)


# Creates Seurat Object ----
# https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
SO <- CreateSeuratObject(counts = count_table, min.cells = 3, min.features  = 100)
SO <- NormalizeData(SO)
SO@meta.data['cell.type'] <- meta_data[,2]
SO@meta.data['time.point'] <- meta_data[,3]
SO@meta.data['tSNE_1'] <- meta_data[,4]
SO@meta.data['tSNE_2'] <- meta_data[,5]
SO@meta.data['nGene'] <- meta_data[,6]
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


# Alternative regression
# Regress out the difference between the G2M and S phase score
# This will maintain the signals separating non-cycling and cycling cells
# but the difference in cell cycle phase of proliferating cells is regressed
# SO$CC.Difference <- SO$S.Score - SO$G2M.Score
# SO <- ScaleData(SO, vars.to.regress = "CC.Difference", features = rownames(SO))

#Save the Obj.
saveRDS(SO, file = sprintf("Scaled_merged_data.rds", dn));
SaveH5Seurat(SO, filename = "Scaled_merged_data.h5Seurat")
Convert("Scaled_merged_data.h5Seurat", dest = "h5ad")



# Differential Expression across time points -----
#Sets cell identity based on amputation
SO <- readRDS('Data/Scaled_merged_data.rds')
Idents(object = SO) <- SO@meta.data$time.point

# SO <- StashIdent(object = SO, save.name = 'Identity')
# SO <- SetAllIdent(object = SO, id = 'orig.ident')

#Finds top DE genes
SO.markers <- FindAllMarkers(object = SO, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
SO.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
top10 <- SO.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)

#Save file
hm <- DoHeatmap(object = SO, features = top10$gene, label = TRUE, size=2,
) + theme(text = element_text(size=4), axis.text.y = element_text(size=4))

jpeg(file="DE_over_time_heatmap.jpeg",
    height = 2800, width = 2300, res=1000);
hm
dev.off()
#CT vs. 2k = sorted vs unsorted









# Looks at DE of oncogenes specifically
onc <- read.csv('Data/Axolotl_oncogenes.csv')
SO.onc <- FindAllMarkers(object = SO, features = onc.genes,
                         only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
SO.onc %>% group_by(cluster) %>% top_n(2, avg_log2FC)
top10onc <- SO.onc %>% group_by(cluster) %>% top_n(10, avg_log2FC)
#NONE OF THE ONCOGENES WERE DIFF EXPRESSED IN THE ENTIRE SUBSET!







# Construct DESeq DataSet  -----
dds <- DESeqDataSetFromMatrix(
  countData = count_table,
  colData = meta_data,
  design = ~orig.ident,
  tidy=FALSE
)
dds

#Sets the level so that the log2 fold changes are calculated
#as amputation over control
dds$orig.ident <- relevel(dds$orig.ident, 'CTa')
dds <- DESeq(dds)

