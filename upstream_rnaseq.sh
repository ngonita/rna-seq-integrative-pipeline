#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# ==============================================================================
# RNA-SEQ UPSTREAM PIPELINE
# FASTQ → BAM → COUNT MATRIX
# ==============================================================================

DIR_DATA="./data/raw_fastq"
DIR_RES="./results"
LOG_DIR="./logs"

GENOME_IDX="./reference/hg38_index"
ANNOTATION="./reference/hg38_annotation.gtf"

THREADS=4

mkdir -p $DIR_RES/qc $DIR_RES/bam $DIR_RES/counts $LOG_DIR

LOG="$LOG_DIR/pipeline.log"
exec > >(tee -a $LOG) 2>&1

echo "Pipeline started at $(date)"

# ---- Tool versions (traceability) ----
fastqc --version
hisat2 --version
samtools --version
featureCounts -v

# ---- Input validation ----
if [ ! -d "$DIR_DATA" ]; then
    echo "ERROR: FASTQ directory not found"
    exit 1
fi

# ---- Sample loop ----
for R1 in $DIR_DATA/*_R1.fastq.gz
do
    SAMPLE=${R1%_R1.fastq.gz}
    BASE=$(basename $SAMPLE)

    echo "Processing sample: $BASE"

    # 1. FASTQC
    fastqc \
        ${SAMPLE}_R1.fastq.gz \
        ${SAMPLE}_R2.fastq.gz \
        -o $DIR_RES/qc \
        -t $THREADS

    # 2. TRIMMOMATIC
    trimmomatic PE -phred33 \
        ${SAMPLE}_R1.fastq.gz \
        ${SAMPLE}_R2.fastq.gz \
        $DIR_RES/${BASE}_R1_clean.fq.gz \
        $DIR_RES/${BASE}_R1_unpaired.fq.gz \
        $DIR_RES/${BASE}_R2_clean.fq.gz \
        $DIR_RES/${BASE}_R2_unpaired.fq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36

    # 3. ALIGNMENT
    hisat2 -p $THREADS \
        -x $GENOME_IDX \
        -1 $DIR_RES/${BASE}_R1_clean.fq.gz \
        -2 $DIR_RES/${BASE}_R2_clean.fq.gz \
        -S $DIR_RES/bam/${BASE}.sam

    # SAM → SORTED BAM
    samtools view -bS $DIR_RES/bam/${BASE}.sam | \
        samtools sort -o $DIR_RES/bam/${BASE}_sorted.bam

    samtools index $DIR_RES/bam/${BASE}_sorted.bam
    rm $DIR_RES/bam/${BASE}.sam

done

# ---- FEATURECOUNTS ----
featureCounts \
    -T $THREADS \
    -p -B -C \
    -t exon \
    -g gene_id \
    -a $ANNOTATION \
    -o $DIR_RES/counts/final_counts_matrix.txt \
    $DIR_RES/bam/*_sorted.bam

echo "Pipeline completed successfully at $(date)"
