#!/bin/bash
# BactPipe - Universal Bacterial Genome Analysis Pipeline
# Usage: bash run.sh

SAMPLE="sample"
THREADS=8

echo "========================================="
echo "BactPipe - Bacterial Analysis Pipeline"
echo "========================================="

# Create directories
mkdir -p data/raw_fastq results/{01_fastqc,02_trimmed,03_assembly,04_quast,05_prokka,06_amr,07_virulence,08_mlst}

# Check if FASTQ files exist
if [ ! -f "data/raw_fastq/${SAMPLE}_R1.fastq" ] && [ ! -f "data/raw_fastq/${SAMPLE}_R1.fastq.gz" ]; then
    echo "ERROR: Please place your FASTQ files in data/raw_fastq/"
    echo "Files should be named: ${SAMPLE}_R1.fastq and ${SAMPLE}_R2.fastq"
    exit 1
fi

echo "[1/8] Quality Control..."
fastqc data/raw_fastq/${SAMPLE}_R1.fastq data/raw_fastq/${SAMPLE}_R2.fastq -o results/01_fastqc -t $THREADS

echo "[2/8] Trimming Reads..."
fastp -i data/raw_fastq/${SAMPLE}_R1.fastq -I data/raw_fastq/${SAMPLE}_R2.fastq \
      -o results/02_trimmed/${SAMPLE}_R1_trimmed.fastq \
      -O results/02_trimmed/${SAMPLE}_R2_trimmed.fastq \
      --detect_adapter_for_pe --cut_front --cut_tail \
      --cut_mean_quality 20 --length_required 50 --thread $THREADS

echo "[3/8] Genome Assembly..."
spades.py -1 results/02_trimmed/${SAMPLE}_R1_trimmed.fastq \
          -2 results/02_trimmed/${SAMPLE}_R2_trimmed.fastq \
          -o results/03_assembly --isolate -t $THREADS -m 32

echo "[4/8] Assembly Quality..."
quast.py results/03_assembly/contigs.fasta -o results/04_quast -t $THREADS

echo "[5/8] Genome Annotation..."
prokka results/03_assembly/contigs.fasta --outdir results/05_prokka \
       --prefix $SAMPLE --kingdom Bacteria --cpus $THREADS --force

echo "[6/8] AMR Detection..."
abricate results/03_assembly/contigs.fasta --db card > results/06_amr/amr_card.tsv

echo "[7/8] Virulence Detection..."
abricate results/03_assembly/contigs.fasta --db vfdb > results/07_virulence/virulence.tsv

echo "[8/8] MLST Typing..."
mlst results/03_assembly/contigs.fasta > results/08_mlst/mlst.txt

echo "========================================="
echo "PIPELINE COMPLETE!"
echo "========================================="
echo "Results: results/assembly/contigs.fasta"
echo "Annotation: results/05_prokka/sample.gbk"
echo "AMR genes: results/06_amr/amr_card.tsv"
echo "MLST: results/08_mlst/mlst.txt"
