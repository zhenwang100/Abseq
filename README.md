# Single-cell 5' RNA sequencing of Bactrian camel PBMCs

#### cellranger

Build reference database and gene expression matrix from raw sequencing data.

#### seurat

Perform QC, integration, dimensional reduction and clustering. 

Use monocle2 to infer differentiation trajectories.

#### trust4

Construct single-cell BCR sequences.

#### Igblast

Annotate BCR sequences and summarize at cell level. 

#### DEG

Differential expression analysis based on pseudobulk method.

Use DESeq2 to identify DEGs and clusterProfiler for function enrichment analysis.

#### annotation_data

Annotation of cell types and BCRs.

- cell_annotation.csv: annotation of cell types for C1-C4 (annotation cluster 9 – T-bet+ B; cluster 10 – naive B; cluster 11 – memory B; cluster 12 – intermediate B)
- cell_immune_annotation.csv: annotation of cell types for C3-C4 during immunization.
- BCR_annotation: annotation of BCR types for C1-C4.
- BCR_immune_annotation.csv: annotation of BCR types for C3-C4 during immunization.
- VHH_sequence: complete VHH sequences assembled for C1-C4.

