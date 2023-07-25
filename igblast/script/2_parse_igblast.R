# Extract V/D/J call, CDR3 and FWR2 from igblast table
args<-commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
	print("Usage: Rscript 2_parse_igblast.R filein fileout")
	q(save = "no")
}
filein<-args[1]		# input file from igblast (AIRR format)
fileout<-args[2]	# output file

# Test file
#filein<-"igblast.txt"
#fileout<-"igblast_parse.txt"

data<-read.table(filein, header = T, sep = "\t")

# Filter sequence without vdj match or chain type
data<-subset(data, sequence_id != '')
data<-subset(data, locus == "VH" | locus == "VL" | locus == "VK")

# vdj start and end
# Note: leading exon is included in v_sequence_start
vdj_start<-data$v_sequence_start + 1	# convert start to 1-based
vdj_end<-data$v_sequence_end

# if d and j genes exist, adjust vdj end
vdj_end[data$d_call != '']<-data$d_sequence_end[data$d_call != '']
vdj_end[data$j_call != '']<-data$j_sequence_end[data$j_call != '']

# Judge VH/VHH based on FR2
# Note: most VH/VHH with complete FR2 [17aa] can be determined with the rule
v_type<-rep(NA, times = nrow(data))
v_type[data$locus == "VH" & nchar(as.character(data$fwr2_aa)) == 17 & grepl("..........GL.....", data$fwr2_aa)]<-"VH"
v_type[data$locus == "VH" & nchar(as.character(data$fwr2_aa)) == 17 & grepl("..........[EQ][R].....", data$fwr2_aa)]<-"VHH"

# Output
# Note: CDR3 annotation is always complete by igblast)
out<-subset(data, select = c("sequence_id", "locus", "productive", "v_call", "d_call", "j_call", "cdr3_aa", "fwr2_aa"))
out<-data.frame(out, v_type = v_type, vdj_start = vdj_start, vdj_end = vdj_end)
write.table(x = out, file = fileout, sep = "\t", row.names = F, quote = F)

