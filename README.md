cat > README.md << 'EOF'
# BactPipe - Bacterial Genome Analysis Pipeline

## One command. Complete analysis. Any bacteria.

### What it does:
- Quality Control (FastQC)
- Read Trimming (fastp)
- Genome Assembly (SPAdes)
- Assembly Quality (QUAST)
- Genome Annotation (Prokka)
- AMR Gene Detection (ABRicate)
- Virulence Factor Detection (ABRicate)
- MLST Typing

### Quick Start

```bash
# 1. Create directory structure
mkdir -p ~/bactpipe_work/Data ~/bactpipe_work/out_results

# 2. Clone pipeline
cd ~/bactpipe_work
git clone https://github.com/Bukhari2917/BactPipe.git

# 3. Place your FASTQ files in Data/ folder
#    Files must be named: sample_R1.fastq and sample_R2.fastq

# 4. Run pipeline
cd BactPipe
bash run.sh
