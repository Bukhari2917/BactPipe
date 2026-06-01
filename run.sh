cd ~/bactpipe_work/BactPipe

cat > run.sh << 'EOF'
#!/bin/bash
# BactPipe - Bacterial Analysis Pipeline
# Usage: bash run.sh

SAMPLE="sample"
THREADS=8

DATA_DIR="../Data"
RESULTS_DIR="../out_results"

echo "========================================="
echo "BactPipe - Bacterial Analysis Pipeline"
echo "========================================="

# Check if FASTQ files exist
if [ ! -f "$DATA_DIR/${SAMPLE}_R1.fastq" ]; then
    echo "ERROR: Cannot find $DATA_DIR/${SAMPLE}_R1.fastq"
    echo ""
    echo "Directory structure:"
    echo "  bactpipe_work/"
    echo "  ├── Data/          ← Place FASTQ files here"
    echo "  │   ├── sample_R1.fastq"
    echo "  │   └── sample_R2.fastq"
    echo "  ├── out_results/   ← Results appear here"
    echo "  └── BactPipe/      ← Pipeline code"
    echo ""
    exit 1
fi

# Create result directories
mkdir -p $RESULTS_DIR/{01_fastqc,02_trimmed,03_assembly,04_quast,05_prokka,06_amr,07_virulence,08_mlst}

echo "[1/8] Quality Control..."
fastqc $DATA_DIR/${SAMPLE}_R1.fastq $DATA_DIR/${SAMPLE}_R2.fastq -o $RESULTS_DIR/01_fastqc -t $THREADS

echo "[2/8] Trimming Reads..."
fastp -i $DATA_DIR/${SAMPLE}_R1.fastq -I $DATA_DIR/${SAMPLE}_R2.fastq \
      -o $RESULTS_DIR/02_trimmed/${SAMPLE}_R1_trimmed.fastq \
      -O $RESULTS_DIR/02_trimmed/${SAMPLE}_R2_trimmed.fastq \
      --detect_adapter_for_pe --cut_front --cut_tail \
      --cut_mean_quality 20 --length_required 50 --thread $THREADS

echo "[3/8] Genome Assembly..."
spades.py -1 $RESULTS_DIR/02_trimmed/${SAMPLE}_R1_trimmed.fastq \
          -2 $RESULTS_DIR/02_trimmed/${SAMPLE}_R2_trimmed.fastq \
          -o $RESULTS_DIR/03_assembly --isolate -t $THREADS -m 32

echo "[4/8] Assembly Quality..."
quast.py $RESULTS_DIR/03_assembly/contigs.fasta -o $RESULTS_DIR/04_quast -t $THREADS

echo "[5/8] Genome Annotation..."
~/prokka-1.14.6/bin/prokka $RESULTS_DIR/03_assembly/contigs.fasta \
       --outdir $RESULTS_DIR/05_prokka \
       --prefix $SAMPLE --kingdom Bacteria --cpus $THREADS --force

echo "[6/8] AMR Detection..."
abricate $RESULTS_DIR/03_assembly/contigs.fasta --db card > $RESULTS_DIR/06_amr/amr_card.tsv

echo "[7/8] Virulence Detection..."
abricate $RESULTS_DIR/03_assembly/contigs.fasta --db vfdb > $RESULTS_DIR/07_virulence/virulence.tsv

echo "[8/8] MLST Typing..."
mlst $RESULTS_DIR/03_assembly/contigs.fasta > $RESULTS_DIR/08_mlst/mlst.txt

echo "========================================="
echo "PIPELINE COMPLETE!"
echo "========================================="
echo "Results saved in: $RESULTS_DIR/"
echo ""
echo "For circular genome visualization:"
echo "  File: $RESULTS_DIR/05_prokka/sample.gbk"
echo "  Upload to: https://proksee.ca/"
echo "========================================="
EOF

chmod +x run.sh
