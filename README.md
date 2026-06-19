# 🧬 BactPipe - Bacterial Whole Genome Analysis Pipeline

**BactPipe** is an automated, end-to-end pipeline for analyzing bacterial whole-genome sequencing (WGS) data. It takes raw sequencing reads (FASTQ files) and produces a comprehensive set of results including genome assembly, gene annotation, antibiotic resistance genes, virulence factors, and MLST typing.

---

## 📖 Table of Contents
1. [What is Whole Genome Sequencing Pipeline?](#-what-is-a-whole-genome-sequencing-pipeline)
2. [What BactPipe Does](#-what-bactpipe-does)
3. [Pipeline Overview](#-pipeline-overview)
4. [Installation](#-installation)
5. [Directory Structure](#-directory-structure)
6. [Preparing Your Data](#-preparing-your-data)
7. [Running the Pipeline](#-running-the-pipeline)
8. [Understanding the Results](#-understanding-the-results)
9. [Troubleshooting](#-troubleshooting)
10. [Citation](#-citation)

---

## 🧬 What is a Whole Genome Sequencing Pipeline?

A **Whole Genome Sequencing (WGS) pipeline** is a series of automated computational steps that process raw DNA sequencing data to extract biological meaning.

### The Journey of Your Data:

```
Raw FASTQ Files (Sequencing Machine Output)
                ↓
         Quality Control (Check read quality)
                ↓
         Read Trimming (Remove bad bases/adapters)
                ↓
         Genome Assembly (Put pieces together like a puzzle)
                ↓
         Assembly Evaluation (Check if genome is complete)
                ↓
         Gene Annotation (Find all genes and their functions)
                ↓
         Downstream Analysis (AMR, Virulence, MLST)
                ↓
         Final Results (Genome, Annotations, Reports)
```

### Why Use a Pipeline?

| Without Pipeline | With Pipeline |
|------------------|---------------|
| Run each tool manually | One command runs everything |
| Hard to reproduce | Fully automated and consistent |
| Takes hours/days | Takes 25-50 minutes |
| Easy to make mistakes | Standardized workflow |
| Hard to track what you did | Logs everything automatically |

---

## 🎯 What BactPipe Does

BactPipe performs **9 key analyses** on your bacterial sequencing data:

| Step | Analysis | Tool | Output |
|------|----------|------|--------|
| 1 | Quality Control | FastQC | HTML quality reports |
| 2 | Read Trimming | Trimmomatic | Cleaned reads (no adapters) |
| 3 | Genome Assembly | SPAdes | Complete genome contigs |
| 4 | Assembly Quality | QUAST | N50, L50, coverage stats |
| 5 | Gene Annotation | Prokka | Annotated genome (.gbk, .gff) |
| 6 | AMR Detection | RGI/CARD | Antibiotic resistance genes |
| 7 | Virulence Detection | Abricate | Virulence factor genes |
| 8 | MLST Typing | mlst | Sequence type (ST) |
| 9 | Functional Annotation | EggNOG-mapper | COG/KOG functions |

**FIX INCLUDED:** This pipeline automatically fixes the common Prokka error *"Prokka needs blastp 2.2 or higher"*.

---

## 🔄 Pipeline Overview

```
                    ┌─────────────────────────────────────────────────────────┐
                    │                                                         │
                    │                  RAW FASTQ FILES                        │
                    │              (sample_R1.fastq, sample_R2.fastq)        │
                    └─────────────────────┬───────────────────────────────────┘
                                          │
                                          ▼
                    ┌─────────────────────────────────────────────────────────┐
                    │  STEP 1: Quality Control (FastQC)                      │
                    │  → Checks read quality scores                          │
                    │  → Identifies adapter contamination                    │
                    │  → Flags overrepresented sequences                     │
                    └─────────────────────┬───────────────────────────────────┘
                                          │
                                          ▼
                    ┌─────────────────────────────────────────────────────────┐
                    │  STEP 2: Read Trimming (Trimmomatic)                   │
                    │  → Removes low-quality bases                           │
                    │  → Removes adapter sequences                           │
                    │  → Filters short reads                                 │
                    └─────────────────────┬───────────────────────────────────┘
                                          │
                                          ▼
                    ┌─────────────────────────────────────────────────────────┐
                    │  STEP 3: Genome Assembly (SPAdes)                      │
                    │  → Assembles reads into contigs                        │
                    │  → Uses k-mer approach                                 │
                    │  → Produces final genome sequence                      │
                    └─────────────────────┬───────────────────────────────────┘
                                          │
                                          ▼
                    ┌─────────────────────────────────────────────────────────┐
                    │  STEP 4: Assembly Quality (QUAST)                      │
                    │  → Calculates N50, L50 values                         │
                    │  → Counts total contigs                                │
                    │  → Checks GC% and coverage                             │
                    └─────────────────────┬───────────────────────────────────┘
                                          │
                                          ▼
                    ┌─────────────────────────────────────────────────────────┐
                    │  STEP 5: Gene Annotation (Prokka)                      │
                    │  → Finds all genes (CDS, rRNA, tRNA)                  │
                    │  → Predicts protein functions                          │
                    │  → Creates .gbk and .gff files                        │
                    └─────────────────────┬───────────────────────────────────┘
                                          │
                                          ▼
        ┌─────────────────────────────────┼─────────────────────────────────┐
        │                                 │                                 │
        ▼                                 ▼                                 ▼
┌───────────────┐               ┌───────────────┐               ┌───────────────┐
│ STEP 6: AMR   │               │ STEP 7:       │               │ STEP 8: MLST  │
│ Detection     │               │ Virulence     │               │ Typing        │
│ (RGI/CARD)    │               │ Detection     │               │ (mlst)        │
│               │               │ (Abricate)    │               │               │
│ Finds genes   │               │ Finds genes   │               │ Finds         │
│ that cause    │               │ that cause    │               │ sequence      │
│ antibiotic    │               │ disease       │               │ type (ST)     │
│ resistance    │               │ (pathogenicity)│              │               │
└───────────────┘               └───────────────┘               └───────────────┘
        │                                 │                                 │
        └─────────────────────────────────┼─────────────────────────────────┘
                                          │
                                          ▼
                    ┌─────────────────────────────────────────────────────────┐
                    │  STEP 9: Functional Annotation (EggNOG-mapper)         │
                    │  → Assigns COG/KOG functional categories               │
                    │  → Identifies biological pathways                      │
                    │  → Provides gene function summaries                    │
                    └─────────────────────┬───────────────────────────────────┘
                                          │
                                          ▼
                    ┌─────────────────────────────────────────────────────────┐
                    │                    FINAL RESULTS                       │
                    │                                                         │
                    │  ✓ Complete genome assembly (contigs.fasta)            │
                    │  ✓ Annotated genome (.gbk, .gff, .faa)                │
                    │  ✓ AMR genes list                                      │
                    │  ✓ Virulence genes list                                │
                    │  ✓ MLST sequence type                                  │
                    │  ✓ Functional annotation (COG/KOG)                     │
                    │  ✓ Quality reports (HTML)                              │
                    └─────────────────────────────────────────────────────────┘
```

---

## 💻 Installation

### Step 1: Install Micromamba

**Why Micromamba?** It's faster, uses less disk space, and works identically to Conda/Mamba.

```bash
# Download and install Micromamba
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Run Micromamba shell setup
./bin/micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc

# Verify installation
micromamba --version
```

![Micromamba Installation](images/micromamba_install.png)

---

### Step 2: Install Git (If Not Installed)

```bash
# For Ubuntu/Debian
sudo apt-get install git

# For CentOS/RHEL
sudo yum install git

# Verify installation
git --version
```

---

### Step 3: Create Directory Structure

```bash
# Create main work directory and subfolders
mkdir -p ~/bactpipe_work/Data
mkdir -p ~/bactpipe_work/out_results

# Go to work directory
cd ~/bactpipe_work
```

---

### Step 4: Download BactPipe

```bash
# Clone the pipeline from GitHub
git clone https://github.com/Bukhari2917/BactPipe.git

# Your final structure will look like:
# ~/bactpipe_work/
# ├── Data/              ← Place your FASTQ files here
# ├── out_results/       ← All results will go here
# └── BactPipe/          ← Pipeline code
#     ├── run.sh
#     ├── environment.yml
#     └── scripts/
```

---

### Step 5: Install Pipeline Dependencies

```bash
# Go to pipeline folder
cd ~/bactpipe_work/BactPipe

# Create conda environment with all required tools
micromamba env create -f environment.yml

# Activate the environment
micromamba activate bactpipe

# Verify all tools are installed correctly
fastqc --version
spades.py --version
prokka --version
abricate --version
rgi --version
mlst --version
```

![Environment Installation](images/env_install.png)

---

### Step 6: Download Databases (One Time)

```bash
# Make sure you're in the bactpipe environment
micromamba activate bactpipe

# Setup Abricate databases (for virulence detection)
abricate --setupdb

# Download EggNOG database (optional - for functional annotation)
# This is large (~10GB) - skip if you don't need EggNOG
eggnog-mapper -m download
```

---

## 📂 Directory Structure Explained

### Complete Folder Structure After Installation:

```
~/bactpipe_work/
│
├── Data/                        # ← PUT YOUR FASTQ FILES HERE
│   ├── sample_R1.fastq          # Forward reads
│   └── sample_R2.fastq          # Reverse reads
│
├── out_results/                 # ← ALL RESULTS WILL APPEAR HERE
│   ├── 01_fastqc/               # Quality control reports
│   ├── 02_trimmed/              # Cleaned reads after trimming
│   ├── 03_assembly/             # Genome assembly files
│   ├── 04_quast/                # Assembly quality metrics
│   ├── 05_prokka/               # Gene annotations
│   ├── 06_amr/                  # AMR genes
│   ├── 07_virulence/            # Virulence factors
│   ├── 08_mlst/                 # MLST typing
│   └── 09_eggnog/               # Functional annotation
│
└── BactPipe/                    # ← PIPELINE CODE (DO NOT MODIFY)
    ├── run.sh                   # Main pipeline script
    ├── environment.yml          # Software dependencies
    ├── scripts/                 # Individual step scripts
    │   ├── 01_fastqc.sh
    │   ├── 02_trim.sh
    │   ├── 03_assembly.sh
    │   ├── 04_quast.sh
    │   ├── 05_prokka.sh
    │   ├── 06_amr.sh
    │   ├── 07_virulence.sh
    │   ├── 08_mlst.sh
    │   └── 09_eggnog.sh
    └── config/                  # Configuration files
        └── adapters.fa          # Adapter sequences for trimming
```

### What Each Directory Contains:

| Directory | Purpose | What You'll Find |
|-----------|---------|------------------|
| `Data/` | Input files | Your raw FASTQ reads |
| `out_results/` | All outputs | Results from every analysis step |
| `BactPipe/` | Pipeline code | Scripts and configuration (don't modify) |

---

## 📂 Preparing Your Data

### File Naming Convention

**CRITICAL:** Your FASTQ files MUST be named exactly as shown below:

| File | Name Format | Example |
|------|-------------|---------|
| Forward reads | `sample_R1.fastq` | `sample_R1.fastq` |
| Reverse reads | `sample_R2.fastq` | `sample_R2.fastq` |

> **Note:** Replace `sample` with your actual sample name (e.g., `Ecoli_R1.fastq`, `S_aureus_R1.fastq`)

### Supported File Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| Uncompressed | `.fastq` | Standard format |
| Uncompressed | `.fq` | Alternative extension |
| Compressed | `.fastq.gz` | Gzipped files (auto-detected) |
| Compressed | `.fq.gz` | Gzipped alternative |

### Place Your Files

```bash
# Copy and rename your files
cp /path/to/your/forward_reads.fastq ~/bactpipe_work/Data/sample_R1.fastq
cp /path/to/your/reverse_reads.fastq ~/bactpipe_work/Data/sample_R2.fastq

# OR if your files are compressed:
cp /path/to/your/forward_reads.fastq.gz ~/bactpipe_work/Data/sample_R1.fastq.gz
cp /path/to/your/reverse_reads.fastq.gz ~/bactpipe_work/Data/sample_R2.fastq.gz

# Verify files are in place
ls -la ~/bactpipe_work/Data/
```

**Expected Output:**
```
-rw-r--r-- 1 user user 398M sample_R1.fastq
-rw-r--r-- 1 user user 398M sample_R2.fastq
```

---

## 🚀 Running the Pipeline

### Step 1: Activate the Environment

```bash
# Always activate this environment before running
micromamba activate bactpipe
```

### Step 2: Run the Pipeline

```bash
# Go to pipeline folder
cd ~/bactpipe_work/BactPipe

# Run the pipeline
bash run.sh
```

### What You Will See

```
=========================================
BactPipe - Bacterial Analysis Pipeline
=========================================
[1/9] Quality Control (FastQC)...     (2-3 minutes)
[2/9] Trimming Reads (Trimmomatic)...  (1-2 minutes)
[3/9] Genome Assembly (SPAdes)...      (15-30 minutes)
[4/9] Assembly Quality (QUAST)...      (1-2 minutes)
[5/9] Genome Annotation (Prokka)...    (2-3 minutes)
[6/9] AMR Detection (RGI)...           (1-2 minutes)
[7/9] Virulence Detection (Abricate).. (1-2 minutes)
[8/9] MLST Typing (mlst)...            (30 seconds)
[9/9] Functional Annotation (EggNOG).. (3-5 minutes)
=========================================
PIPELINE COMPLETE!
=========================================
Total time: 25-50 minutes (depends on genome size)
```

### Runtime Estimates

| Genome Size | Time | RAM Usage |
|-------------|------|-----------|
| Small (~2 Mbp) | 15-25 min | 4-8 GB |
| Medium (~4 Mbp) | 25-40 min | 8-12 GB |
| Large (~6 Mbp) | 40-60 min | 12-16 GB |

---

## 📊 Understanding the Results

All results are saved in `~/bactpipe_work/out_results/`:

### Complete Results Structure

```
out_results/
│
├── 01_fastqc/                          # Quality Control Reports
│   ├── sample_R1_fastqc.html          # ← Open in browser
│   ├── sample_R1_fastqc.zip
│   ├── sample_R2_fastqc.html
│   └── sample_R2_fastqc.zip
│
├── 02_trimmed/                         # Cleaned Reads
│   ├── sample_R1_paired.fastq         # High-quality forward reads
│   ├── sample_R1_unpaired.fastq       # Orphan forward reads
│   ├── sample_R2_paired.fastq         # High-quality reverse reads
│   └── sample_R2_unpaired.fastq       # Orphan reverse reads
│
├── 03_assembly/                        # Genome Assembly
│   ├── contigs.fasta                  # ← YOUR GENOME SEQUENCE
│   ├── scaffolds.fasta
│   └── assembly_graph/
│
├── 04_quast/                           # Assembly Quality
│   ├── report.html                    # ← Open in browser
│   ├── report.txt
│   ├── contigs_report.tsv
│   └── coverage/
│
├── 05_prokka/                          # Genome Annotation
│   ├── sample.gbk                     # ← Annotated genome (for visualization)
│   ├── sample.gff                     # Gene features (for IGV)
│   ├── sample.faa                     # Protein sequences
│   ├── sample.fna                     # Nucleotide sequences
│   ├── sample.tsv                     # Gene table (Excel-friendly)
│   ├── sample.sqn                     # For NCBI submission
│   └── sample.txt                     # Summary statistics
│
├── 06_amr/                             # AMR Genes
│   └── amr_card.tsv                   # ← Antibiotic resistance genes
│
├── 07_virulence/                       # Virulence Factors
│   └── virulence.tsv                  # ← Virulence genes found
│
├── 08_mlst/                            # MLST Typing
│   └── mlst.txt                       # ← Sequence type (ST)
│
└── 09_eggnog/                          # Functional Annotation
    ├── eggnog_results.tsv             # ← COG/KOG functions
    └── eggnog_summary.txt             # Summary statistics
```

---

### Key Result Files and How to Use Them

| File | What It Contains | How to Open/Use |
|------|------------------|-----------------|
| **contigs.fasta** | Your complete genome sequence | Open in any text editor, or upload to NCBI BLAST |
| **report.html** | Assembly quality metrics (N50, L50, # contigs) | Double-click to open in any web browser |
| **sample.gbk** | Annotated genome with all genes and features | Open in Geneious, SnapGene, or Proksee for visualization |
| **sample.gff** | Gene coordinates and functional annotations | Import into IGV (Integrative Genomics Viewer) |
| **sample.faa** | All predicted protein sequences | Use for BLAST searches against other proteins |
| **amr_card.tsv** | List of antibiotic resistance genes found | Open in Excel or text editor |
| **virulence.tsv** | List of virulence factor genes found | Open in Excel or text editor |
| **mlst.txt** | Sequence type (ST) number | Compare with outbreak strains |
| **eggnog_results.tsv** | COG/KOG functional categories of all genes | Open in Excel or text editor |

---

### Quick View Commands

```bash
# View AMR genes
cat ~/bactpipe_work/out_results/06_amr/amr_card.tsv

# View MLST type
cat ~/bactpipe_work/out_results/08_mlst/mlst.txt

# View virulence genes
cat ~/bactpipe_work/out_results/07_virulence/virulence.tsv

# View assembly stats
cat ~/bactpipe_work/out_results/04_quast/report.txt

# Count how many contigs
grep -c "^>" ~/bactpipe_work/out_results/03_assembly/contigs.fasta
```

---

### How to Check If Your Assembly is Good

Open `~/bactpipe_work/out_results/04_quast/report.html` in your browser and look for:

| Metric | Good | Excellent | Explanation |
|--------|------|-----------|-------------|
| **N50** | > 50,000 bp | > 100,000 bp | Half the genome is in contigs this long |
| **# Contigs** | < 100 | < 50 | Fewer contigs = more complete genome |
| **GC%** | Matches species | Matches species | Should be typical for your organism |
| **Total Length** | Matches expected | Matches expected | Should match known genome size |

---

## 🛠️ Troubleshooting

### Common Issues and Solutions

| Issue | Solution |
|-------|----------|
| `micromamba: command not found` | Run `source ~/.bashrc` or restart terminal |
| `environment.yml not found` | Make sure you're in the BactPipe folder: `cd ~/bactpipe_work/BactPipe` |
| `Prokka needs blastp 2.2 or higher` | **Already fixed in this pipeline!** The fix is included. |
| `SPAdes fails with out of memory` | SPAdes needs ~8-16 GB RAM. Close other programs or use a smaller dataset. |
| `Pipeline stops at EggNOG` | Check internet connection. EggNOG downloads databases on first run. |
| `RGI not found` | Run `micromamba install -c conda-forge rgi` |
| `Abricate database missing` | Run `abricate --setupdb` |
| `Permission denied` | Make scripts executable: `chmod +x BactPipe/scripts/*.sh` |

### Fix Common Errors

```bash
# If micromamba environment fails to activate
micromamba init bash
source ~/.bashrc

# If a tool is missing
micromamba install -c conda-forge -c bioconda [tool_name]

# To completely reset the environment
micromamba env remove -n bactpipe
micromamba env create -f environment.yml

# To rerun just one step (example: run Prokka again)
cd ~/bactpipe_work/BactPipe
bash scripts/05_prokka.sh
```

### Check Logs

```bash
# View pipeline logs
tail -f ~/bactpipe_work/BactPipe/logs/*.log
```

---

## 📝 Citation

If you use BactPipe in your research, please cite:

```
Bukhari, S.A. (2026). BactPipe: An automated bacterial whole-genome analysis pipeline.
GitHub Repository: https://github.com/Bukhari2917/BactPipe
```

### BibTeX Entry

```bibtex
@misc{bactpipe2026,
  author = {Bukhari, S.A.},
  title = {BactPipe: An automated bacterial whole-genome analysis pipeline},
  year = {2026},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/Bukhari2917/BactPipe}}
}
```

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

BactPipe relies on these wonderful open-source tools:

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - Quality control
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) - Read trimming
- [SPAdes](https://cab.spbu.ru/software/spades/) - Genome assembly
- [QUAST](http://quast.sourceforge.net/) - Assembly evaluation
- [Prokka](https://github.com/tseemann/prokka) - Genome annotation
- [RGI / CARD](https://card.mcmaster.ca/) - AMR detection
- [Abricate](https://github.com/tseemann/abricate) - Virulence detection
- [mlst](https://github.com/tseemann/mlst) - MLST typing
- [EggNOG-mapper](http://eggnog-mapper.embl.de/) - Functional annotation

---

**Happy Analyzing! 🧬
