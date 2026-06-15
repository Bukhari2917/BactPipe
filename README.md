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

### Step 1: Install Micromamba

```bash
# Download and install Micromamba
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Run Micromamba shell setup
./bin/micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc
