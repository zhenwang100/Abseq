# Run monocle2 pipeline (R4.0 only supports monocle2)
library(Seurat)
library(monocle)

# Only extract B cell clusters from Seurat object
dir<-Sys.getenv(x = "PROJ")
load(dir, "/seurat/integrated/comb.Rdata")
B_data<-subset(comb_data, idents = c('0', '1', '5', '8', '13'))

# Construct input data for monocle2
# Expression matrix (RNA count is recommended)
expr<-B_data$RNA@counts

# Cell annotation (using Seurat clusters)
p_data<-B_data@meta.data 
p_data$celltype<-B_data@active.ident  

# Gene annotation (gene_short_name is required)
f_data<-data.frame(gene_short_name = row.names(expr), row.names = row.names(expr))

pd<-new('AnnotatedDataFrame', data = p_data) 
fd<-new('AnnotatedDataFrame', data = f_data)
cds<-newCellDataSet(expr, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

# Estimate size factor and dispersion
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)

head(pData(cds))		# size factor for cells
head(dispersionTable(cds))	# dispersion for genes

# Filter low-expressed genes
cds<-detectGenes(cds, min_expr = 0.1)
expressed_genes<-row.names(subset(fData(cds), num_cells_expressed >= 10))

head(fData(cds))	# num_cells_expressed for genes

# Identify DEGs for ordering cells
diff<-differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 20)
diff_genes<-subset(diff, qval < 0.01)
diff_genes<-diff_genes[order(diff_genes$qval, decreasing = F),]

order_genes<-rownames(diff_genes)[1:3000]
cds<-setOrderingFilter(cds, order_genes)
ggsave("order_gene_DEG.pdf", plot_ordering_genes(cds))

# Order cells
cds<-reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds<-orderCells(cds)

head(pData(cds))	# Pseudotime and state for cells
table(pData(cds)$celltype, pData(cds)$State)	# Summary cell type by state

# Plot trajectory by cell type
traj_celltype_DEG<-plot_cell_trajectory(cds, color_by = "celltype", show_backbone = TRUE)
ggsave("traj_celltype_DEG.pdf", traj_celltype_DEG)

# Plot trajectory by state
traj_state_DEG<-plot_cell_trajectory(cds, color_by = "State", show_backbone = TRUE)
ggsave("traj_state_DEG.pdf", traj_state_DEG)

# Plot trajectory by pseudotime (root state should be assigned according to data)
cds<-orderCells(cds, root_state = 4)
traj_pseudo_DEG<-plot_cell_trajectory(cds, color_by = "Pseudotime", show_backbone = TRUE)
ggsave("traj_pseudo_DEG.pdf", traj_pseudo_DEG)

# Save monocle dataset
save(cds, file = "cds.RData")

# Perform branch-specific gene selection analysis (branch_point should be assigned according to data)
BEAM_res<-BEAM(cds[order_genes,], branch_point = 2)
BEAM_res<-BEAM_res[order(BEAM_res$qval),]
BEAM_res<-BEAM_res[,c("gene_short_name", "pval", "qval")]

pdf("BEAM_heatmap.pdf")
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res, qval < 1e-4)),], branch_point = 2, num_clusters = 6, cores = 1, use_gene_short_name = T, show_rownames = T)
dev.off()

# Plot marker genes
markers<-c("CR2", "ITGAX", "TBX21", "CXCR4", "CD69", "SELL", "ITGB1", "ITGB7", "S100A6", "S100A7", "S100A13", "IGHA", "IGHM", "CD38")
pdf("marker_heatmap.pdf")
plot_genes_branched_heatmap(cds[markers,], branch_point = 2, num_clusters = 4, cores = 1, use_gene_short_name = T, show_rownames = T)
dev.off()

