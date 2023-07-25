# Surport old OS
export TENX_IGNORE_DEPRECATED_OS=1

# Mapping
cd $PROJ/cellranger

cellranger count --id=C1 --transcriptome=ref --fastqs=$PROJ/data/C1 --expect-cells 5000
cellranger count --id=C2 --transcriptome=ref --fastqs=$PROJ/data/C2 --expect-cells 5000

cellranger count --id=C3 --transcriptome=ref --fastqs=$PROJ/data/C3 --expect-cells 5000
cellranger count --id=C4 --transcriptome=ref --fastqs=$PROJ/data/C4 --expect-cells 5000

cellranger count --id=C3_42d --transcriptome=ref --fastqs=$PROJ/data/C3_42d --expect-cells 5000
cellranger count --id=C3_42d --transcriptome=ref --fastqs=$PROJ/data/C3_42d --expect-cells 5000

cellranger count --id=C4_56d --transcriptome=ref --fastqs=$PROJ/data/C4_56d --expect-cells 5000
cellranger count --id=C4_56d --transcriptome=ref --fastqs=$PROJ/data/C4_56d --expect-cells 5000

