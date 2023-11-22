setwd("/projects/academic/rgunawan/Panos/10X")
dirSel = "/projects/academic/rgunawan/Panos/x86_64-pc-linux-gnu-library/4.0"

# Add the path to your library search path
.libPaths( c(dirSel ,.libPaths()))

# load libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(scales)
library(readxl)

# files paths
decont_file='/projects/academic/rgunawan/Panos/10X/data_tissue/decontXcounts_df_glp.csv'
log_norm_file='/projects/academic/rgunawan/Panos/10X/data_tissue/X_unscaled_matrix_df_glp.csv'
meta_matrix_file='/projects/academic/rgunawan/Panos/10X/data_tissue/Meta_matrix_glp.csv'
cell_ids_file='/projects/academic/rgunawan/Panos/10X/data_tissue/cell_names_all_glp.csv'
glyco_groups='/projects/academic/rgunawan/Panos/10X/GlycoEnzOntoShort.xlsx'

# load data 
# mat_raw is the raw counts, mat_log is the log1p normalized data, and mat_meta has the metadata dataframes
mat_raw <- read.csv(decont_file,header = TRUE)
mat_meta <- read.csv(meta_matrix_file,header = TRUE)
mat_log<- read.csv(log_norm_file,header = TRUE)
ro_nam<- read.csv(cell_ids_file,header = FALSE)
ronam1=as.vector(t(ro_nam))

# Set the barcodes as the row names
rownames(mat_raw) <- ronam1
rownames(mat_meta) <- ronam1
rownames(mat_log) <- ronam1

# Transpose the data and convert to sparse matrix
mat_raw <- as(t(as.matrix(mat_raw)), "sparseMatrix")
mat_log <- as(t(as.matrix(mat_log)), "sparseMatrix")

# Create the Seurat object
seurat_1 <- CreateSeuratObject(counts = mat_raw, meta.data = mat_meta)
seurat_1@assays$RNA@data <- mat_log

# Define the groups that are going to be compared for the DE analysis
# The following DE analysis is for glycopathways / tissues.
# In Tabula Sapiens, value = "tissue_in_publication" for tissues, and value = "cell_types" for cell types.
seurat_1 <- SetIdent(seurat_1, value = "tissue_in_publication")

# Perform DE analysis. 
# Output: p-values, avgLFC, clusters (tissues) and glycopathways that meet
# the criteria of logfc.threshold (logfc threshold) and return.thresh (p-values threshold)
tissue_markers_all_pv01 <- FindAllMarkers(
  seurat_1,
  #only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.00001,
  test.use ='MAST',
  return.thresh=0.1
)

# Save the result
# markers_file1='/projects/academic/rgunawan/Panos/10X/data_tissue_1/markers_all_tissues_mast_pv01_nolfc.csv'
# write.csv(tissue_markers_all_pv01,markers_file1)