# Run in conda env R4.0
source ~/.bashrc
conda activate R4.0

if [ $# -ne 1 ]; then
	echo "Usage: sh run_trust4.sh SAMPLE"
	exit 1
fi

#SAMPLE=C1
SAMPLE=$1	

run-trust4 -f $PROJ/trust4/ref/IR_ref.fa --ref $PROJ/trust4/ref/IR_ref.fa \
 -1 $PROJ/data/$SAMPLE/*_R1_*.fastq.gz \
 -2 $PROJ/data/$SAMPLE/*_R2_*.fastq.gz \
 -o $SAMPLE \
 --od $PROJ/trust4/$SAMPLE \
 --barcode $PROJ/data/$SAMPLE/*_R1_*.fastq.gz \
 --barcodeRange 0 15 + \
 --UMI $PROJ/data/$SAMPLE/*_R1_*.fastq.gz \
 --umiRange 16 25 + \
 --read1Range 40 -1 \
 -t 30 \

# Exception when add barcode list for unknown reasons... 
#--barcodeWhitelist $HOME/prog/install/cellranger-6.1.2/lib/python/cellranger/barcodes/737K-august-2016.txt

