# Run on R4.0 and Seurat4.0 
# Examine one sample
library(Seurat)
library(ggplot2)

# Usage: Rscript one_sample.R C1
args<-commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
	print("Usage: Rscript one_sample.R sample")
	q(save = "no")
}

# Get project directory
dir<-Sys.getenv(x = "PROJ")

sample<-args[1]
data_dir<-paste(dir, "/cellranger/", sample, "/outs/filtered_feature_bc_matrix", sep = "")

# Output directory
out_dir<-paste(dir, "/seurat/", sample, sep = "")
setwd(out_dir)

# Read data
sample_data<-Read10X(data.dir = data_dir)
sample_data<-CreateSeuratObject(counts = sample_data, project = sample, assay = "RNA")

# MT gene content
sample_data[["percent.mt"]]<-PercentageFeatureSet(sample_data, 
	features = c('ND1', 'ND2', 'COX1', 'COX2', 'ATP8', 'ATP6', 'COX3', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'CYTB'))

# SCT normalization
sample_data<-SCTransform(sample_data, verbose = FALSE)

# content of SCT assay
#sample_data[["SCT"]]@counts		# corrected counts
#sample_data[["SCT"]]@data		# log1p(corrected counts)
#sample_data[["SCT"]]@scale.data	# pearson residuals
#sample_data[["SCT"]]@meta.features	# Feature-level meta data (eg. sct mean, variance of each gene)
#sample_data[["SCT"]]@var.features	# Selected variable features (3000 by default)

# Dimensinal reduction
sample_data<-RunPCA(sample_data, verbose = FALSE)
sample_data<-RunUMAP(sample_data, dims = 1:30, verbose = FALSE)

# Clustering
sample_data<-FindNeighbors(sample_data, dims = 1:30, verbose = FALSE)
sample_data<-FindClusters(sample_data, verbose = FALSE)

# Plot cluster
cluster<-DimPlot(sample_data, label = TRUE)
ggsave("cluster.pdf", cluster)

# Plot UMI count and mt gene content across cluster
sample_data[["log_nCount_RNA"]]<-log10(sample_data[["nCount_RNA"]])
qc<-FeaturePlot(sample_data, features = c("log_nCount_RNA", "percent.mt"))
ggsave("qc.pdf", qc, width = 14, height = 7)

# Plot markers
markers<-DotPlot(sample_data, features = c(
	"CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "LOC102522090",	# T cells (LOC102522090, TRDC)
	"KLRB1", "NKG7", "LOC102516303", 		# NK cells (LOC102516303, CD16)
	"CD19", "MS4A1", "CD38",	# B cells
	"CD14", "CD68",			# Monocytes
	"ITGAX", "FCER1A",		# cDC
	"IL3RA", "NRP1"			# pDC
	),) + RotatedAxis()	
ggsave("markers.pdf", markers)



