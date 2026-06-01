# BactPipe - Bacterial Genome Analysis Pipeline

**One command. Complete analysis. Circular genome visualization.**

## Quick Start

```bash
# 1. Create directories
mkdir -p ~/bactpipe_work/Data ~/bactpipe_work/out_results

# 2. Clone pipeline
cd ~/bactpipe_work
git clone https://github.com/Bukhari2917/BactPipe.git

# 3. Place your FASTQ files (rename to sample_R1.fastq, sample_R2.fastq)
cp /path/to/your_R1.fastq Data/sample_R1.fastq
cp /path/to/your_R2.fastq Data/sample_R2.fastq

# 4. Run pipeline (ONE COMMAND)
cd BactPipe
bash run.sh
