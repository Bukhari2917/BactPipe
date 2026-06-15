#!/bin/bash

# BactPipe - Bacterial Analysis Pipeline
# Author: Bukhari2917

set -e

echo "========================================="
echo "BactPipe - Bacterial Analysis Pipeline"
echo "========================================="

# Check if Data directory has FASTQ files
if [ ! -d "../Data" ]; then
    echo "ERROR: Data directory not found!"
    exit 1
fi

# Create results directory
mkdir -p ../out_results

# Step 1: Quality Control
echo "[1/8] Quality Control..."
fastqc ../Data/*.fastq -o ../out_results/01_fastqc/ -t 8

# Step 2: Trimming
echo "[2/8] Trimming Reads..."
fastp -i ../Data/sample_R1.fastq -I ../Data/sample_R2.fastq \
      -o ../out_results/02_trimmed/trimmed_R1.fastq.gz \
      -O ../out_results/02_trimmed/trimmed_R2.fastq.gz \
      --html ../out_results/02_trimmed/fastp_report.html \
      --thread 8

# Step 3: Assembly
echo "[3/8] Genome Assembly..."
spades.py -1 ../out_results/02_trimmed/trimmed_R1.fastq.gz \
          -2 ../out_results/02_trimmed/trimmed_R2.fastq.gz \
          -o ../out_results/03_assembly/ \
          --isolate --threads 8

# Step 4: Fix Prokka (ADD THIS SECTION)
echo "[4/8] Fixing Prokka for compatibility..."
if [ -f "scripts/fix_prokka.sh" ]; then
    source scripts/fix_prokka.sh
else
    # Direct fix if script not found
    PROKKA_PATH=$(which prokka)
    sed -i 's/MINVER  => "2.2"/MINVER  => "2.0"/g' "$PROKKA_PATH"
fi

# Step 5: Annotation
echo "[5/8] Genome Annotation..."
prokka --outdir ../out_results/05_prokka/ \
       --prefix sample \
       --kingdom Bacteria \
       --addgenes \
       --cpus 8 \
       --force \
       ../out_results/03_assembly/contigs.fasta

# Step 6: AMR Detection
echo "[6/8] AMR Detection..."
abricate ../out_results/03_assembly/contigs.fasta > ../out_results/06_amr/amr_card.tsv

# Step 7: Virulence Detection
echo "[7/8] Virulence Detection..."
abricate --db vfdb ../out_results/03_assembly/contigs.fasta > ../out_results/07_virulence/virulence.tsv

# Step 8: MLST Typing
echo "[8/8] MLST Typing..."
mlst ../out_results/03_assembly/contigs.fasta > ../out_results/08_mlst/mlst.txt

echo "========================================="
echo "PIPELINE COMPLETE!"
echo "Results are in: ../out_results/"
echo "========================================="
