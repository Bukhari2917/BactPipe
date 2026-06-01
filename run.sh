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
    echo "Your directory structure should be:"
    echo "  bactpipe_work/"
    echo "  ├── Data/          ← Place your FASTQ files here"
    echo "  │   ├── sample_R1.fastq"
    echo "  │   └── sample_R2.fastq"
    echo "  ├── out_results/   ← Results will appear here"
    echo "  └── BactPipe/      ← Pipeline code (you are here)"
    echo ""
    exit 1
fi

mkdir -p $RESULTS_DIR/{01_fastqc,02_trimmed,03_assembly,04_quast,05_prokka,06_amr,07_virulence,08_mlst,09_visualizations}

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

# Fix GBK for circular genome
echo "========================================="
echo "FIXING GBK FOR CIRCULAR GENOME"
echo "========================================="

CLEAN_DIR="$RESULTS_DIR/05_prokka_clean"
FINAL_GBK="$RESULTS_DIR/05_prokka/sample_fixed.gbk"

mkdir -p $CLEAN_DIR

echo "[1/3] Creating clean contig names..."
sed 's/ .*//; s/_cov_[0-9.]*//g; s/_length_[0-9]*//g' $RESULTS_DIR/03_assembly/contigs.fasta > $CLEAN_DIR/contigs_clean.fasta

echo "[2/3] Re-annotating with Prokka..."
~/prokka-1.14.6/bin/prokka $CLEAN_DIR/contigs_clean.fasta \
       --outdir $CLEAN_DIR \
       --prefix sample \
       --kingdom Bacteria \
       --cpus $THREADS \
       --force \
       --locustag SAMPLE

echo "[3/3] Copying GBK file..."
cp $CLEAN_DIR/sample.gbk $FINAL_GBK

echo ""
echo "========================================="
echo "CIRCULAR GENOME READY!"
echo "========================================="
echo "Fixed GBK file: $FINAL_GBK"
echo ""
echo "Download to your computer:"
echo "  scp sayed@j3-053-010:$FINAL_GBK ./"
echo ""
echo "Then upload to: https://proksee.ca/"
echo "========================================="
