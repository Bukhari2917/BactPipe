# 🧬 BactPipe - Bacterial Analysis Pipeline

**BactPipe** is an automated pipeline that analyzes bacterial DNA sequencing data. It takes your raw FASTQ files and produces:

- Complete bacterial genome assembly
- Annotated genome with all genes
- Antibiotic resistance genes (AMR)
- Virulence factors
- Bacterial strain type (MLST)
- Circular genome visualization file

---

## 📋 Table of Contents
- [Installation](#-installation)
- [Preparing Your Data](#-preparing-your-data)
- [Running the Pipeline](#-running-the-pipeline)
- [Results](#-results)
- [Troubleshooting](#-troubleshooting)

---

## 💻 Installation

### Step 1: Install Micromamba

```bash
# Download and install Micromamba
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Run Micromamba shell setup
./bin/micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc
```

![Micromamba Installation](images/micromamba_install.png)

### Step 2: Create Directory Structure

```bash
# Create main work directory and subfolders
mkdir -p ~/bactpipe_work/Data
mkdir -p ~/bactpipe_work/out_results

# Go to work directory
cd ~/bactpipe_work
```

### Step 3: Download BactPipe

```bash
# Clone the pipeline
git clone https://github.com/Bukhari2917/BactPipe.git

# Your structure will look like:
# bactpipe_work/
# ├── Data/
# ├── out_results/
# └── BactPipe/
```

### Step 4: Install Dependencies with Micromamba

```bash
# Go to pipeline folder
cd BactPipe

# Create conda environment using the environment.yml file
micromamba env create -f environment.yml

# Activate the environment
micromamba activate bactpipe

# Verify installation
fastqc --version
spades.py --version
prokka --version
```

---
## 📂 Preparing Your Data

### File Naming Convention

**IMPORTANT**: Your FASTQ files MUST be named exactly as:

| File | Name Format |
|------|-------------|
| Forward reads | `sample_R1.fastq` |
| Reverse reads | `sample_R2.fastq` |

### Place Your Files

```bash
# Copy and rename your files
cp /path/to/your/forward_reads.fastq ~/bactpipe_work/Data/sample_R1.fastq
cp /path/to/your/reverse_reads.fastq ~/bactpipe_work/Data/sample_R2.fastq

# Verify files are in place
ls -la ~/bactpipe_work/Data/
```

---

## 🚀 Running the Pipeline

```bash
# Go to pipeline folder
cd ~/bactpipe_work/BactPipe

# Activate the environment
micromamba activate bactpipe

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
Total time: 25-50 minutes
```

---

## 📊 Results

All results are saved in `~/bactpipe_work/out_results/`:

```
out_results/
├── 01_fastqc/          # Quality control reports
├── 02_trimmed/         # Cleaned reads
├── 03_assembly/        # Genome assembly (contigs.fasta)
├── 04_quast/           # Assembly quality report
├── 05_prokka/          # Genome annotation (.gbk, .gff, .faa)
├── 06_amr/             # AMR genes
├── 07_virulence/       # Virulence factors
├── 08_mlst/            # MLST typing
└── 09_eggnog/          # Functional annotation
```

### Key Files

| File | What It Contains |
|------|------------------|
| `contigs.fasta` | Your assembled genome |
| `report.html` | Assembly quality (N50, L50) |
| `sample.gbk` | Annotated genome for visualization |
| `amr_card.tsv` | Antibiotic resistance genes |
| `virulence.tsv` | Virulence genes |
| `mlst.txt` | Sequence type (ST) |

---

## 🛠️ Troubleshooting

| Issue | Solution |
|-------|----------|
| `micromamba: command not found` | Run `source ~/.bashrc` |
| `Prokka needs blastp 2.2 or higher` | Already fixed in this pipeline! |
| Out of memory | SPAdes needs ~8-16 GB RAM |

---

## 📝 Citation

If you use BactPipe, please cite:

```
Bukhari, S.A. (2026). BactPipe: An automated bacterial whole-genome analysis pipeline.
GitHub: https://github.com/Bukhari2917/BactPipe
```

---

## 📜 License

MIT License - see [LICENSE](LICENSE) file for details.

---

Happy Analyzing!
