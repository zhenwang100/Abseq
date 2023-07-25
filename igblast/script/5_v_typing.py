# Summarize v types (including VH/VHH, CDR3 and VL/VK) by barcodes
# Note c type is only included for barcodes with VH/VHH or CDR3 annotations
import sys

if len(sys.argv) != 3:
	sys.exit("Usage: python 5_v_typing.py filein fileout") 

filein = sys.argv[1]	# Input file from 4_merge_vdjc.R
fileout = sys.argv[2]	# Output file

# Test file
#filein = "vdjc_merge.txt"
#fileout = "v_type.txt"

barcodes = set()	# Unique barcode set

contig_locus = {}	# Map from contig to locus (VH/VL/VK)
contig_fwr2 = {}	# Map from contig to FWR2
contig_cdr3 = {}	# Map from contig to CDR3
contig_c = {}		# Map from contig to c type

barcode_vh_contig = {}		# Map from barcode to VH contig (if multiple contigs, choose the one with the highest coverage)
barcode_vhh_contig = {}		# Map from barcode to VHH contig (see above)
barcode_cdr3_contig = {}	# Map from barcode to CDR3 contig (see above)
barcode_c_contig = {}		# Map from barcode to c gene contig (see above)
barcode_l_contig = {}		# Map from barcode to L/K contig (see above

barcode_vh_cov = {}		# Map from barcode to highest VH contig coverage
barcode_vhh_cov = {}		# Map from barcode to highest VHH contig coverage
barcode_cdr3_cov = {}		# Map from barcode to highest CDR3 contig coverage
barcode_c_cov = {}		# Map from barcode to highest c gene contig coverage
barcode_l_cov = {}		# Map from barcode to highest L/K contig coverage

filein = open(filein)
fileout = open(fileout, "w")
fileout.write("barcode\tvh_contig\tvhh_contig\tvh_cov\tvhh_cov\tv_type\tfwr2\tcdr3_contig\tcdr3\tc_contig\tc_type\tl_contig\tl_type\tl_cdr3\n")

filein.readline()	# head line
for line in filein:
	cols = line.strip().split("\t")
	contig = cols[0]
	barcode = contig.split("_")[0]
	cov = float(cols[2])
	locus = cols[3]
	cdr3 = cols[8]
	fwr2 = cols[9]
	vtype = cols[10]
	ctype = cols[13]

	# Assign contig-level map
	contig_locus[contig] = locus
	contig_fwr2[contig] = fwr2
	contig_cdr3[contig] = cdr3
	contig_c[contig] = ctype

	# Assign barcode-level map
	if locus == 'VH':
		# Select VH contig with the highest coverage
		if vtype == 'VH' and (barcode not in barcode_vh_contig or barcode_vh_cov[barcode] < cov):
			barcode_vh_contig[barcode] = contig
			barcode_vh_cov[barcode] = cov
			barcodes.add(barcode)
		# Select VHH contig with the highest coverage
		if vtype == 'VHH' and (barcode not in barcode_vhh_contig or barcode_vhh_cov[barcode] < cov):
			barcode_vhh_contig[barcode] = contig
			barcode_vhh_cov[barcode] = cov
			barcodes.add(barcode)
		# Select CDR3 contig with the highest coverage
		if cdr3 != '' and (barcode not in barcode_cdr3_contig or barcode_cdr3_cov[barcode] < cov):
			barcode_cdr3_contig[barcode] = contig
			barcode_cdr3_cov[barcode] = cov
			barcodes.add(barcode)
		# Select C gene contig with the highest coverage
		if ctype != 'NA' and (barcode not in barcode_c_contig or barcode_c_cov[barcode] < cov):
			barcode_c_contig[barcode] = contig
			barcode_c_cov[barcode] = cov
			barcodes.add(barcode)
	# Select K/L contig with the highest coverage
	if locus == "VK" or locus == "VL":
		if barcode not in barcode_l_contig or barcode_l_cov[barcode] < cov:
			barcode_l_contig[barcode] = contig
			barcode_l_cov[barcode] = cov

# Match barcode-level information
for barcode in barcodes:
	(vh_contig, vhh_contig, vh_cov, vhh_cov, v_type, fwr2, cdr3_contig, cdr3, c_contig, c_type, l_contig, l_type, l_cdr3) = ('NA',) * 13
	# VH/VHH contigs and coverage
	if barcode in barcode_vh_contig:
		vh_contig = barcode_vh_contig[barcode]
		vh_cov = barcode_vh_cov[barcode]
	if barcode in barcode_vhh_contig:
		vhh_contig = barcode_vhh_contig[barcode]
		vhh_cov = barcode_vhh_cov[barcode]

	# Assign v type by comparing VH and VHH coverage
	if (vhh_contig == 'NA' and vh_contig != 'NA') or (vhh_contig != 'NA' and vh_contig != 'NA' and vh_cov > vhh_cov):
		v_type = 'VH'
		fwr2 = contig_fwr2[vh_contig]
	if (vh_contig == 'NA' and vhh_contig != 'NA') or (vh_contig != 'NA' and vhh_contig != 'NA' and vhh_cov > vh_cov):
		v_type = 'VHH'
		fwr2 = contig_fwr2[vhh_contig]

	# CDR3 contigs and sequences
	if barcode in barcode_cdr3_contig:
		cdr3_contig = barcode_cdr3_contig[barcode]
		cdr3 = contig_cdr3[cdr3_contig]

	# C gene contigs and sequences
	if barcode in barcode_c_contig:
		c_contig = barcode_c_contig[barcode]
		c_type = contig_c[c_contig]

	# L/K contigs and CDR3 sequences
	if barcode in barcode_l_contig:
		l_contig = barcode_l_contig[barcode]
		l_type = contig_locus[l_contig]
		l_cdr3 = contig_cdr3[l_contig]
		
	# Output
	fileout.write("\t".join((barcode, vh_contig, vhh_contig, str(vh_cov), str(vhh_cov), v_type, fwr2, cdr3_contig, cdr3, c_contig, c_type, l_contig, l_type, l_cdr3)) + "\n")

