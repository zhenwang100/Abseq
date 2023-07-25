# Get possible junction C genes
import re
import sys

if len(sys.argv) != 3:
	sys.exit("Usage: python 3_parse_cblast.py filein fileout") 

filein = sys.argv[1]	# input blast file (table format)
fileout = sys.argv[2]	# output file

# Test file
#filein = "cblast.txt"
#fileout = "cblast_parse.txt"

filein = open(filein)
fileout = open(fileout, "w")
fileout.write("sequence_id\tc_call\tc_sequence_start\tc_sequence_end\tc_germline_start\tc_germline_end\n")

# Current query sequence
cur_seqid = None

for line in filein:
	cols = line.strip().split("\t")
	seqid = cols[0]		# query sequence
	c_call = cols[1]        # target sequence
	c_seq_start = cols[6]   # query start
	c_seq_end = cols[7]     # query end
	c_ref_start = cols[8]   # target start
	c_ref_end = cols[9]     # target end

	# Adjust C gene names with highly similar sequences
	if re.search("/[MS]", c_call):          # neglect M/S type for IGHC 
		c_call = c_call.split("/")[0]
	if c_call == "IGHG1B1":
		c_call = "IGHG1B"
	if re.match("IGLC", c_call):            # neglect IGLC gene number 
		c_call = "IGLC"

	# Consider best hit only
	if cur_seqid is None or cur_seqid != seqid:
		cur_seqid = seqid	# Refresh
		# Write junction C genes (match from first 3 bp in refrence)
		# Note: <= 3 (allow for blast mismatch at beginning) and == 1 (most strict) do not diff too much
		if int(c_ref_start) <= 3:
			fileout.write("\t".join((seqid, c_call, c_seq_start, c_seq_end, c_ref_start, c_ref_end)) + "\n")

