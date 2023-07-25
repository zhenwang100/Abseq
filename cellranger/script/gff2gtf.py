# This script convert NCBI gff to gtf for cellranger.
# Unlike other public tools like gffread, or gt gff3_to_gtf, this script use ID directly from original gff.
import sys
import gzip
import re

# Input: NCBI gff file, gzip
# Output: gtf file
if len(sys.argv) != 3:
	print("Usage: python gff2gtf.py gff_file.gz gtf_file")
	sys.exit()

gff_file = sys.argv[1]
gtf_file = sys.argv[2]

gff = gzip.open(gff_file)
gtf = open(gtf_file, 'w')

# Current gene and transcript
gene_id = None
gene_name = None
gene_biotype = None
transcript_id = None
exon_id = None

# Summary all possible biotypes
gene_biotypes = {}

# For each line
for line in gff:
	if re.match("#", line):
		continue
	cols = line.strip().split("\t")

	# Gene entry
	if cols[2] == 'gene':
		attr = cols[8]
		# gene_id
		m = re.match("ID=([^;]+)", attr)
		if m is not None:
			gene_id = m.group(1)
		else:
			print("Error: cannot find ID at line:" + line)
			sys.exit()
		# gene_name
		m = re.search("Name=([^;]+)", attr)
		if m is not None:
			gene_name = m.group(1)
		else:
			print("Error: cannnot find Name at line:" + line)
			sys.exit()
		# gene_biotype
		m = re.search("gene_biotype=([^;]+)", attr)
		if m is not None:
			gene_biotype = m.group(1)
			# Summary biotypes
			if gene_biotype not in gene_biotypes:
				gene_biotypes[gene_biotype] = 0
			else:
				gene_biotypes[gene_biotype] += 1
		else:
			print("Error: cannot find gene_biotype at line:" + line)
			sys.exit()
		# gtf gene entry
		new_attr = "gene_id \"" + gene_id + "\"; gene_name \"" + gene_name + "\"; gene_biotype \"" + gene_biotype + "\";"
		cols[8] = new_attr
		gtf.write("\t".join(cols) + "\n")

		
	# Transcript entry (many types of transcript...)
	if cols[2] in ('mRNA', 'transcript', 'C_gene_segment', 'V_gene_segment', 'lncRNA', 'rRNA', 'tRNA', 'snRNA', 'snoRNA'):
		attr = cols[8]

		# Check gene parent
		m = re.search("Parent=([^;]+)", attr)
		if m is None or m.group(1) != gene_id:
			continue

		# transcript_id
		m = re.match("ID=([^;]+)", attr)
		if m is not None:
			transcript_id = m.group(1)
		else:
			print("Error: cannot find ID at line:" + line)
			sys.exit()
		# gtf transcript entry
		new_attr = "gene_id \"" + gene_id + "\"; gene_name \"" + gene_name + "\"; gene_biotype \"" + gene_biotype + "\"; transcript_id \"" + transcript_id + "\";"
		cols[2] = 'transcript'
		cols[8] = new_attr
		gtf.write("\t".join(cols) + "\n")

	if cols[2] == 'exon':
		attr = cols[8]

		# Check transcript parent
		m = re.search("Parent=([^;]+)", attr)
                if m is None or m.group(1) != transcript_id:
                        continue

		# exon_id
		m = re.match("ID=([^;]+)", attr)
                if m is not None:
                        exon_id = m.group(1)
                else:
                        print("Error: cannot find ID at line:" + line)
                        sys.exit()
		# gtf exon entry
		new_attr = "gene_id \"" + gene_id + "\"; gene_name \"" + gene_name + "\"; gene_biotype \"" + gene_biotype + "\"; transcript_id \"" + transcript_id + "\"; exon_id \"" + exon_id + "\";"
		cols[8] = new_attr
		gtf.write("\t".join(cols) + "\n")

print("Successfully finish. Gene_biotype summary:")
print(gene_biotypes)

