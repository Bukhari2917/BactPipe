# BactPipe - Bacterial Analysis Pipeline

BactPipe is an automated pipeline that analyzes bacterial DNA sequencing data. It takes your raw FASTQ files and produces:
- Complete bacterial genome assembly
- Annotated genome with all genes
- Antibiotic resistance genes (AMR)
- Virulence factors
- Bacterial strain type (MLST)
- Circular genome visualization file

> **✅ FIX INCLUDED:** This pipeline automatically fixes the common Prokka error *"Prokka needs blastp 2.2 or higher"*.

---

## PART 1: INSTALLATION (One Time Only)

### Step 1: Install Miniconda

```bash
# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install it
bash Miniconda3-latest-Linux-x86_64.sh

# Follow prompts:
# - Press Enter to read license
# - Type 'yes' to accept
# - Press Enter for default location
# - Type 'yes' to initialize

# Reload terminal
source ~/.bashrc
