# Combine v types across samples
#samples<-c("20211102-1", "20211102-2", "14d-5", "14d-6")
samples<-c("20211102-1", "20211102-2", "14d-5", "14d-6", "42d-5", "42d-6", "56d-5", "56d-6")

data<-data.frame()
for (i in 1 : length(samples)) {
	sample<-samples[i]
	sample_file<-paste("/picb/bigdata/project/AbSeq/igblast/", sample, "/v_type.txt", sep = "")
	sample_data<-read.table(sample_file, sep = "\t", header = T)

	# Convert barcode to cellranger ("-1") and then Seurat format ("_i") (i is sample index)
	sample_data$barcode<-paste(sample_data$barcode, 1, sep = "-")
	sample_data$barcode<-paste(sample_data$barcode, i, sep = "_")

	# Combine
	data<-rbind(data, sample_data)
}

write.table(file = "comb_v_type.txt", x = data, sep = "\t", quote = F, row.names = F)

