# R4.0 and Seurat 4.0
# Integrate C3, C3_42d, C3_56d, C4, C4_42d, C4_56d and perform reference-based transfer annotation
library(Seurat)
library(ggplot2)

# Samples and UMI cutoffs
samples<-c("C3", "C3_42d", "C3_56d", "C4", "C4_42d", "C4_56d")
low_cutoff<-c(1000, 1000, 1000, 1000, 1000, 1000)
high_cutoff<-c(60000, 60000, 60000, 60000, 60000, 60000)

names(low_cutoff)<-samples
names(high_cutoff)<-samples

# Output directory
dir<-Sys.getenv(x = "PROJ")
out_dir<-paste(dir, "/seurat/integrated_immune/")
setwd(out_dir)

# Read data and data QC
data_list<-lapply(samples, FUN = function(sample) {
	# read data
	data_dir<-paste(dir, "/cellranger/", sample, "/outs/filtered_feature_bc_matrix", sep = "")
	sample_data<-Read10X(data.dir = data_dir)

	# project should be different for each sample
	sample_data<-CreateSeuratObject(counts = sample_data, project = sample, assay = "RNA")

	# data QC
	sample_data[["percent.mt"]]<-PercentageFeatureSet(sample_data, 
		features = c('ND1', 'ND2', 'COX1', 'COX2', 'ATP8', 'ATP6', 'COX3', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'CYTB'))
	sample_data<-SubsetData(sample_data, subset.name = "nCount_RNA", 
		low.threshold = low_cutoff[sample], high.threshold = high_cutoff[sample])
	sample_data<-SubsetData(sample_data, subset.name = "percent.mt", high.threshold = 5)
})

# SCT normalization
data_list<-lapply(data_list, FUN = SCTransform)

# Integrate data based on SCT
features<-SelectIntegrationFeatures(data_list, nfeatures = 3000)
data_list<-PrepSCTIntegration(data_list, anchor.features = features)
anchors<-FindIntegrationAnchors(data_list, normalization.method = "SCT", anchor.features = features)
comb_immune_data<-IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Dimensional reduction
comb_immune_data<-RunPCA(comb_immune_data)
comb_immune_data<-RunUMAP(comb_immune_data, reduction = "pca", dims = 1:30)

# Clustering (default resolution 0.8)
comb_immune_data<-FindNeighbors(comb_immune_data, reduction = "pca", dims = 1:30)
comb_immune_data<-FindClusters(comb_immune_data)

cluster<-DimPlot(comb_immune_data, label = TRUE)
ggsave("cluster_0.8.pdf", cluster)

# Save integrated data
save(list = ls(all.names = TRUE), file = "comb_immune.Rdata")

# PBMC markers
markers<-DotPlot(comb_immune_data, features = c(
	"CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "LOC102522090",	# T cells (LOC102522090: TRDC)
	"KLRB1", "NKG7", "LOC102516303", 		# NK cells (LOC102516303: CD16)
	"CD19", "MS4A1", "CD38",	# B cells
	"CD14", "CD68",			# Monocytes
	"TYMS",				# Profiling cells
	"ITGAX", "FCER1A",		# cDC
	"IL3RA", "NRP1"			# pDC
	),) + RotatedAxis()
ggsave("markers_0.8.pdf", markers)

# Perform reference-based annotation (annotated comb_data)
load("../integrated/comb.Rdata")

DefaultAssay(comb_data)<-"SCT"
comb_data<-FindVariableFeatures(comb_data)

# Query
DefaultAssay(comb_immune_data)<-"SCT"

# Transfer cell cluster labels
anchors<-FindTransferAnchors(reference = comb_data, query = comb_immune_data, normalization.method = "SCT")
ann_cluster<-TransferData(anchorset = anchors, refdata = comb_data$annotation_cluster)
comb_immune_data$annotation_cluster<-ann_cluster$predicted.id

cell_type<-TransferData(anchorset = anchors, refdata = comb_data$cell_type)
comb_immune_data$cell_type<-cell_type$predicted.id

# Plot transferred cell cluster in UMAP of query data
ann_cluster<-DimPlot(comb_immune_data, group.by = 'annotation_cluster', label = TRUE)
ggsave("ann_cluster.pdf", ann_cluster)

cell_type<-DimPlot(comb_immune_data, group.by = 'cell_type', label = TRUE)
ggsave("cell_type.pdf", cell_type)

# Summary count
table(comb_immune_data$orig.ident, comb_immune_data$cell_type)
table(comb_immune_data$orig.ident, comb_immune_data$annotation_cluster)

# Preserve data
save(comb_immune_data, file = "comb_immune.Rdata")

