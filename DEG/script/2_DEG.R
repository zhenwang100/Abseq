library(DESeq2)

# For each cell type or subtype
for (cell_type in c("B", "CD4-T", "CD8-T", "DC", "gdT", "Mono", "NK", "plasma-B", "ProT")) {
	cnt.file<-paste("pseudobulk/", cell_type, "_cnt.txt", sep = "")		# Input file:pseudobulk count
	cpm.file<-paste("pseudobulk/", cell_type, "_cpm.txt", sep = "")		# Output file: CPM normalization
	stat.file<-paste("pseudobulk/", cell_type, "_stat.txt", sep = "")	# Output file: test statistics
	
	# Count data and meta info
	# 'Replicate' is biological replicate, 'label' is before or after immunization
	data<-read.table(cnt.file, header = T)
	meta<-data.frame(replicate = c("C3", "C4", "C3", "C4", "C3", "C4"), 
		label = c("before", "before", "after", "after", "after", "after"))

	dds<-DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~replicate + label)
	
	# CPM normalization
	cpm<-fpm(dds)
	write.table(cpm, file = cpm.file, sep = "\t", quote = F)

	# Pre-filtering low expression genes (only minimum filtering is needed here)
	keep<-rowSums(counts(dds)) >= 10
	dds<-dds[keep,]

	# DEG
	dds<-DESeq(dds)
	res<-results(dds)
	
	res<-as.data.frame(res)
	res<-res[order(res$padj),]

	# padj < 0.05
	write.table(subset(res, padj <= 0.05,), file = stat.file, sep = "\t", quote = F)
	# export all
	#write.table(res, file = stat.file, sep = "\t", quote = F)
}

