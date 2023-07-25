# R4.0 and Seurat 4.0
# Integrate samples C1-C4
library(Seurat)
library(ggplot2)

# Samples and UMI cutoffs
samples<-c("C1", "C2", "C3", "C4")
low_cutoff<-c(1000, 1000, 1000, 1000)
high_cutoff<-c(60000, 60000, 60000, 60000)

names(low_cutoff)<-samples
names(high_cutoff)<-samples

# Output directory
dir<-Sys.getenv(x = "PROJ")
out_dir<-paste(dir, "/seurat/integrated/")
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
comb_data<-IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Dimensional reduction
comb_data<-RunPCA(comb_data)
comb_data<-RunUMAP(comb_data, reduction = "pca", dims = 1:30)

# Clustering (default resolution 0.8)
comb_data<-FindNeighbors(comb_data, reduction = "pca", dims = 1:30)
comb_data<-FindClusters(comb_data)

cluster<-DimPlot(comb_data, label = TRUE)
ggsave("cluster_0.8.pdf", cluster)

# Save integrated data
save(list = ls(all.names = TRUE), file = "comb.Rdata")

# Optional: Change resulution
for (res in c(0.1, 0.3, 0.5, 1.0, 1.2, 1.5)) {
	comb_data<-FindClusters(comb_data, resolution = res)
	pdf<-paste("cluster_", res, ".pdf", sep = "")
	ggsave(pdf, DimPlot(comb_data, label = T))
}
# Choose a resolution
Idents(comb_data)<-comb_data@meta.data$integrated_snn_res.0.8

# Gene level analysis should be performed for SCT assay (not integrated or RNA), including expression visualization
DefaultAssay(comb_data)<-"SCT"

# Plot known markers (using both DotPlot and FeaturePlot is useful)
# PBMC markers
markers<-DotPlot(comb_data, features = c(
	"CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "LOC102522090",	# T cells (LOC102522090: TRDC)
	"KLRB1", "NKG7", "LOC102516303", 		# NK cells (LOC102516303: CD16)
	"CD19", "MS4A1", "CD38",	# B cells
	"CD14", "CD68",			# Monocytes
	"TYMS",				# Profiling cells
	"ITGAX", "FCER1A",		# cDC
	"IL3RA", "NRP1"			# pDC
	),) + RotatedAxis()
ggsave("markers_0.8.pdf", markers)

# Identifying markers for a cluster
markers_0<-FindMarkers(comb_data, ident.1 = 0, ident.2 = c(1, 5, 8), assay = "SCT")

