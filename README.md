cat > README.md << 'EOF'
# BactPipe - Universal Bacterial Genome Analysis Pipeline

**One command. Any bacteria. 8 complete analyses.**

[![Nextflow](https://img.shields.io/badge/nextflow-22.10.6-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## What It Does

| # | Analysis | Tool | Output |
|---|----------|------|--------|
| 1 | Quality Control | FastQC | HTML reports |
| 2 | Read Trimming | fastp | Cleaned FASTQ |
| 3 | Genome Assembly | SPAdes | Contigs/Scaffolds |
| 4 | Assembly Quality | QUAST | Metrics report |
| 5 | Genome Annotation | Prokka | GBK, GFF, FAA |
| 6 | AMR Detection | ABRicate | Resistance genes |
| 7 | Virulence Detection | ABRicate | Virulence factors |
| 8 | MLST Typing | mlst | Sequence type |

## Quick Start

```bash
# 1. Clone repository
git clone https://github.com/Bukhari2917/BactPipe.git
cd BactPipe

# 2. Install dependencies
conda env create -f environment.yml
conda activate bactpipe

# 3. Place your FASTQ files
mkdir -p data/raw_fastq
cp /path/to/your/sample_R1.fastq.gz data/raw_fastq/
cp /path/to/your/sample_R2.fastq.gz data/raw_fastq/

# 4. Run pipeline (ONE COMMAND!)
nextflow run main.nf
