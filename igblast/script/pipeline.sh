# Re-annotate Trust4 contigs using igblast
# Set sample name
if [ $# -ne 1 ]; then
        echo "Usage: sh pipeline.sh SAMPLE"
        exit 1
fi

#SAMPLE=C1
SAMPLE=$1

# Create working directory
WD=$PROJ/igblast/$SAMPLE
if [ ! -d $WD ]; then
 	mkdir $WD
fi

# Reference database (blast database for V/D/J and C)
DB=$PROJ/igblast/ref
# Trust4 contigs
IN=$PROJ/trust4/$SAMPLE/${SAMPLE}_annot.fa
# Igblast output
OUT=$WD/igblast.txt
# Additional output for blast C genes
OUT_C=$WD/cblast.txt

# Perform igblast for vdj annotation
# -auxiliary_data should be provided for CDR3 annotation
igblastn -germline_db_V $DB/IG_V -germline_db_D $DB/IG_D -germline_db_J $DB/IG_J \
	-auxiliary_data $DB/camel_gl.aux \
	-query $IN \
	-ig_seqtype Ig \
	-show_translation -outfmt 19 -out $OUT \
	-num_threads 20 -evalue 1e-5 \
	-num_alignments_V 1

# Perform blast for c annotation
blastall -p blastn -d $DB/IG_C -i $IN -e 1e-5 -m 8 -a 20 -o $OUT_C

# In-house scripts to summarize the results
# Get contig length and coverage
python 1_parse_contig.py $IN $WD/contig.txt

# Get V/D/J call, alignment length, CDR3 and FWR2
Rscript 2_parse_igblast.R $WD/igblast.txt $WD/igblast_parse.txt

# Get C call for junctions within contigs
python 3_parse_cblast.py $WD/cblast.txt $WD/cblast_parse.txt

# Merge V/D/J call and C call and preserve valide contigs only
Rscript 4_merge_vdjc.R $WD/contig.txt $WD/igblast_parse.txt $WD/cblast_parse.txt $WD/vdjc_merge.txt

# Summarize v types by barcodes
python 5_v_typing.py $WD/vdjc_merge.txt $WD/v_type.txt

