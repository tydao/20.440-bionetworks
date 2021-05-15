
# Changes directory ----
dn = "/Users/Dao/Documents/MIT Grad/Classes/20.440/Project"
setwd(dn)

# Loads in the Data----
library(Seurat)
library(SeuratDisk)
library(magrittr)
library(dplyr)
library(ggplot2)


## Function to perform preprocessing on the counts table ----
# Takes in a dataframe of counts and metadata and name of output file
# https://satijalab.org/seurat/articles/essential_commands.html
preprocess_counts_table <- function(df, output_fn) {
  #Creates meta data table
  meta_data <- df[1:6]
  colnames(meta_data) <- colnames(df)[1:6]
  rownames(meta_data) <- meta_data$cell_id
  
  #Creates counts table (gene x library/cell)
  count_table <- df[, !names(df) %in% 
                      c('cell_id','ident', 'orig.ident',
                        'tSNE_1', 'tSNE_2', 'nGene')]
  rownames(count_table) <- df[,1]
  count_table <- t(count_table)
  
  # Creates Seurat Object and filters out cells with less than 100 genes
  # https://scrnaseq-course.cog.sanger.ac.uk/website/seurat-chapter.html
  SO <- CreateSeuratObject(counts = count_table, min.cells = 3, min.features  = 100)
  SO <- NormalizeData(SO, normalization.method = "LogNormalize")
  
  # Adds in meta data
  SO@meta.data['cell.type'] <- meta_data[colnames(SO),2]
  SO@meta.data['time.point'] <- meta_data[colnames(SO),3]
  SO@meta.data['tSNE_1'] <- meta_data[colnames(SO),4]
  SO@meta.data['tSNE_2'] <- meta_data[colnames(SO),5]
  SO@meta.data['nGene'] <- meta_data[colnames(SO),6]
  
  # Subset on the expression level of PRRX1
  #     SO <- subset(x = SO, subset = PRRX1 > 0)
  
  # Find variable features
  SO <- FindVariableFeatures(object = SO, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
  
  # Removes cell-cycle variation ----
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We then
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # Gets cell cycle score and changes the current identity
  SO <- CellCycleScoring(SO, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  RidgePlot(SO, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
  
  #     # "Before" PCA
  #     SO <- RunPCA(SO, features = c(s.genes, g2m.genes))
  #     DimPlot(SO)
  
  # Regressed PCA
  SO <- ScaleData(SO, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SO))
  SO <- RunPCA(SO, features = c(s.genes, g2m.genes))
  #     DimPlot(SO)
  
  #Save the Obj.
  if (!missing(output_fn)){
    saveRDS(SO, file = sprintf(paste(output_fn,".rds", sep="")))
    SaveH5Seurat(SO, filename = paste(output_fn,".h5Seurat", sep=""))
    Convert(paste(output_fn,".h5Seurat", sep=""), dest = "h5ad")
  }

  return (SO)
}

# Run script ----

# Process table 5 (healthy samples)
t5 <- read.csv("Data/aaq0681_TableS5.csv") #Read in csv
# preprocess_counts_table(t5, "Table5_processed")
t5_SO <-preprocess_counts_table(t5)

# Process table 9 (amputations)
t9 <- read.csv('aa10681_TableS9.csv')

t9_18 <- fn[fn$orig.ident == '18dpaA' | fn$orig.ident == '18dpaB',]
preprocess_counts_table(t9_18, "Table9_18dpa_processed")

t9_25 <- fn[fn$orig.ident == '25dpaA' | fn$orig.ident == '25dpaB',]
preprocess_counts_table(t9_25, "Table9_25dpa_processed")

t9_36 <- fn[fn$orig.ident == '36dpaA' | fn$orig.ident == '36dpaB',]
preprocess_counts_table(t9_18, "Table9_36dpa_processed")


