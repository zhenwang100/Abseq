cd $PROJ/cellranger

# Convert NCBI gff to gtf
python script/gff2gtf.py NCBI/GCF_009834535.1_BCGSAC_Cfer_1.0_genomic.gff.gz gtf/tmp.gtf

# Combine mt transcript annotation and manual C gene annotation
cat gtf/tmp.gtf gtf/mt_transcript.gtf gtf/C_gene_segment.gtf > gtf/genome.gtf

# Preserve only protein_coding gene and C gene
cellranger mkgtf gtf/genome.gtf gtf/genome_filter.gtf --attribute=gene_biotype:protein_coding --attribute=gene_biotype:BCR/TCR_gene_segment

# Make reference
cellranger mkref --genome=ref --fasta=NCBI/GCF_009834535.1_BCGSAC_Cfer_1.0_genomic.fna --genes=gtf/genome_filter.gtf

