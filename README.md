# BactPipe - Bacterial Genome Analysis Pipeline

Complete pipeline for bacterial genome analysis including assembly, annotation, AMR, virulence, CRISPR, and MLST.

## Quick Start

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Clone pipeline
git clone https://github.com/Bukhari2917/BactPipe.git
cd BactPipe

# Place your FASTQ files
# Files must be named: sample_R1.fastq and sample_R2.fastq

# Run pipeline
nextflow run main.nf --reads 'data/*_R{1,2}.fastq' -profile conda
