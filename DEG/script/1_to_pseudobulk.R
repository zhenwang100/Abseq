# Use Libra to convert cell matrix to pseudobulk matrix
# However, Libra may have bugs in pseudobulk DEG analysis...
library(stringr)
library(Libra)
library(Seurat)

load("comb_immune.Rdata")

# Sample name
samples<-str_split_fixed(comb_immune_data$orig.ident, "-", 2) %>% data.frame()

# 'label' is experiment condition, 'replicate' is biological replicate
comb_immune_data$label<-samples$X1
comb_immune_data$replicate<-samples$X2

# Use counts from RNA assay (raw), but not SCT assay (corrected)
DefaultAssay(comb_immune_data)<-"RNA"

# 1. Convert to pseudobulks for major cell types
exp<-to_pseudobulk(comb_immune_data)
sapply(names(exp), FUN = function(cell_type) {
	out.file = paste("pseudobulk/", cell_type, "_cnt.txt", sep = '')
	write.table(x = exp[[cell_type]], file = out.file, sep = "\t", quote = F)
})

# Run DE with edgeR (default) or DESeq, but there are bugs in both methods...
#DE<-run_de(comb_immune_data)
#DE<-run_de(comb_immune_data, de_family = 'pseudobulk', de_method = 'DESeq2')
