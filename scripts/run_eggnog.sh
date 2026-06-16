#!/bin/bash
# eggNOG-mapper for functional annotation (COG, KEGG, GO)

set -e

echo "========================================="
echo "Running eggNOG-mapper"
echo "========================================="

# Check if protein file exists
if [ ! -f "../out_results/05_prokka/sample.faa" ]; then
    echo "ERROR: Protein file not found. Run main pipeline first."
    exit 1
fi

# Create output directory
mkdir -p ../out_results/09_eggnog

# Setup database directory
EGGNOG_DB="$HOME/eggnog_db"
mkdir -p $EGGNOG_DB

# Download database if not exists (first time only)
if [ ! -f "$EGGNOG_DB/eggnog.db" ]; then
    echo "Downloading eggNOG database (first time, ~40GB)..."
    cd $EGGNOG_DB
    wget -c http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz && gunzip eggnog.db.gz
    wget -c http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz && gunzip eggnog_proteins.dmnd.gz
    wget -c http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz && tar -xzf eggnog.taxa.tar.gz
fi

# Create environment if not exists
if ! micromamba env list | grep -q "eggnog_env"; then
    micromamba create -n eggnog_env python=3.9 pip -y
    micromamba run -n eggnog_env pip install eggnog-mapper
fi

# Run annotation
cd ../out_results/05_prokka/
micromamba run -n eggnog_env emapper.py -i sample.faa \
    -o ../09_eggnog/functional_annotation \
    --data_dir $EGGNOG_DB \
    --cpu 8 \
    --tax_scope bacteria \
    --override

echo "========================================="
echo "Complete! Results: ../out_results/09_eggnog/"
echo "========================================="
