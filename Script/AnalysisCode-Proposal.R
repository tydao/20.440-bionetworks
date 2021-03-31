library(Seurat)

# Read in the data
raw_counts <- read.csv('Data/aaq0681_TableS5.csv')
#colnames(raw_counts) <- raw_counts[1, ]
#raw_counts <- raw_counts[-1, ]
rownames(raw_counts) <- raw_counts[, 1]
raw_counts <- raw_counts[, -1]
raw_counts[c('ident', 'orig.ident', 'tSNE_1', 'tSNE_2', 'nGene')] <- NULL
SO <- CreateSeuratObject(t(raw_counts))

# Prepare Seurat Data
SO <- NormalizeData(SO)
all.genes <- rownames(SO)
SO <- FindVariableFeatures(object = SO)
SO <- ScaleData(SO, features = all.genes)
SO <- RunPCA(SO, features = VariableFeatures(object = SO))

# Run Clustering and Reduction
SO <- RunUMAP(SO, dims = 1:10)
SO <- FindClusters(SO, reduction.type = "umap", resolution = 0.051)

# Plot Graph
DimPlot(SO, reduction = "umap", label = TRUE)
FeaturePlot(SO, features=c("TNMD", "ASPN", "SPARC", "HMGN2"))
