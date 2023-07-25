# Merge igblast (vdj) and junction c genes
args<-commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
	print("Usage: Rscript 4_merge_vdjc.R contig igblast_parse cblast_parse fileout")
	q(save = "no")
}

contig<-args[1]			# Input file 1 from 1_parse_contig.py
igblast_parse<-args[2]		# Input file 2 from 2_parse_igblast.R
cblast_parse<-args[3]		# Input file 3 from 3_parse_cblast.py
fileout<-args[4]		# Output file

# Test file
#contig<-"contig.txt"
#igblast_parse<-"igblast_parse.txt"
#cblast_parse<-"cblast_parse.txt"
#fileout<-"vdjc_merge.txt"

cblast<-read.table(cblast_parse, sep = "\t", header = T)
igblast<-read.table(igblast_parse, sep = "\t", header = T)
contig<-read.table(contig, sep = "\t", header = T)

# Merge the three tables (left join for igblast)
data<-merge(igblast, cblast, by = "sequence_id", all.x = T)
data<-merge(contig, data, by = "sequence_id")

# Length of vdj and c region
vdj_len<-data$vdj_end - data$vdj_start + 1
c_len<-data$c_sequence_end - data$c_sequence_start + 1

# Index for valide vdj_c junction
jun_index<-!is.na(data$c_call) & data$j_call != '' & abs(data$c_sequence_start - data$vdj_end) <= 3
# Preserve c call for valide junctions only
data$c_call[!jun_index]<-NA
vdj_len[jun_index]<-vdj_len[jun_index] + c_len[jun_index]

# Index for productive vdj
prod_index<-data$productive | is.na(data$productive)

# Length QC (loose cutoff: 150, strict cutff: 200)
data_qc<-subset(data, prod_index & vdj_len > 150)

write.table(x = data_qc, file = fileout, sep = "\t", quote = F, row.names = F)

