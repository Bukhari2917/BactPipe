# BactPipe

[![Nextflow](https://img.shields.io/badge/nextflow-24.04.2-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**BactPipe** is a complete bacterial genome analysis pipeline.

## Quick Start

```bash
nextflow run Bukhari2917/BactPipe --reads 'data/*_R{1,2}.fastq' -profile conda