# Extact contig length and coverage
# Note: the contig length and log(coverage) are highly correlated.
import re
import sys

if len(sys.argv) != 3:
	sys.exit("Usage: python 1_parse_contig.py filein fileout") 

filein = sys.argv[1]	# input contig file (trust4 fasta format)
fileout = sys.argv[2]	# output file

filein = open(filein)
fileout = open(fileout, "w")

fileout.write("sequence_id\tlength\tcoverage\n")
for line in filein:
	m = re.match(">(.+)", line)
	if m is not None:
		title = m.group(1)
		(seqid, length, coverage) = title.split(" ")[0:3]
		fileout.write(seqid + "\t" + length + "\t" + coverage + "\n")

