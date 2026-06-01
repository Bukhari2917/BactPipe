cd ~/bactpipe_work/BactPipe

cat > README.md << 'EOF'
# BactPipe - Bacterial Genome Analysis Pipeline

**One command. Complete analysis.**

## Quick Start

```bash
# 1. Create directories
mkdir -p ~/bactpipe_work/Data ~/bactpipe_work/out_results

# 2. Clone pipeline
cd ~/bactpipe_work
git clone https://github.com/Bukhari2917/BactPipe.git

# 3. Place FASTQ files in Data/ (rename to sample_R1.fastq, sample_R2.fastq)

# 4. Run pipeline
cd BactPipe
bash run.sh
