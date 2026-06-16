#!/bin/bash

# BactPipe - Bacterial Analysis Pipeline
# Author: Bukhari2917

set -e

echo "========================================="
echo "BactPipe - Bacterial Analysis Pipeline"
echo "========================================="

# Check if data directory has FASTQ files
if [ ! -d "../data" ] && [ ! -d "../Data" ]; then
    echo "ERROR: data directory not found!"
    exit 1
fi

# Use data or Data folder
if [ -d "../data" ]; then
    DATA_DIR="../data"
else
    DATA_DIR="../Data"
fi

# Create results directory
mkdir -p ../out_results

# Step 1: Quality Control
echo "[1/9] Quality Control..."
mkdir -p ../out_results/01_fastqc
fastqc $DATA_DIR/*.fastq -o ../out_results/01_fastqc/ -t 8

# Step 2: Trimming
echo "[2/9] Trimming Reads..."
mkdir -p ../out_results/02_trimmed
fastp -i $DATA_DIR/sample_R1.fastq -I $DATA_DIR/sample_R2.fastq \
      -o ../out_results/02_trimmed/trimmed_R1.fastq.gz \
      -O ../out_results/02_trimmed/trimmed_R2.fastq.gz \
      --html ../out_results/02_trimmed/fastp_report.html \
      --thread 8

# Step 3: Assembly
echo "[3/9] Genome Assembly..."
mkdir -p ../out_results/03_assembly
spades.py -1 ../out_results/02_trimmed/trimmed_R1.fastq.gz \
          -2 ../out_results/02_trimmed/trimmed_R2.fastq.gz \
          -o ../out_results/03_assembly/ \
          --isolate --threads 8

# Step 4: Fix Prokka
echo "[4/9] Fixing Prokka for compatibility..."
if [ -f "scripts/fix_prokka.sh" ]; then
    bash scripts/fix_prokka.sh
else
    PROKKA_PATH=$(which prokka)
    sed -i 's/MINVER  => "2.2"/MINVER  => "2.0"/g' "$PROKKA_PATH"
fi

# Step 5: Annotation
echo "[5/9] Genome Annotation..."
mkdir -p ../out_results/05_prokka
prokka --outdir ../out_results/05_prokka/ \
       --prefix sample \
       --kingdom Bacteria \
       --addgenes \
       --cpus 8 \
       --force \
       ../out_results/03_assembly/contigs.fasta

# Step 6: AMR Detection
echo "[6/9] AMR Detection..."
mkdir -p ../out_results/06_amr
abricate ../out_results/03_assembly/contigs.fasta > ../out_results/06_amr/amr_card.tsv
abricate --summary ../out_results/06_amr/amr_card.tsv > ../out_results/06_amr/amr_summary.txt

# Step 7: Virulence Detection
echo "[7/9] Virulence Detection..."
mkdir -p ../out_results/07_virulence
abricate --db vfdb ../out_results/03_assembly/contigs.fasta > ../out_results/07_virulence/virulence.tsv

# Step 8: MLST Typing
echo "[8/9] MLST Typing..."
mkdir -p ../out_results/08_mlst
mlst ../out_results/03_assembly/contigs.fasta > ../out_results/08_mlst/mlst.txt

# Step 9: eggNOG functional annotation
echo "[9/9] Functional Annotation (COG/KEGG/GO)..."
if [ -f "scripts/run_eggnog.sh" ]; then
    bash scripts/run_eggnog.sh
else
    echo "WARNING: scripts/run_eggnog.sh not found, skipping eggNOG"
fi

echo "========================================="
echo "PIPELINE COMPLETE!"
echo "Results are in: ../out_results/"
echo "========================================="
