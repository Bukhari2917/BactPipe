#!/bin/bash
# Fix Prokka BLAST version detection bug
# This script patches Prokka to recognize BLAST 2.16+ as valid

echo "========================================="
echo "Applying Prokka BLAST version fix..."
echo "========================================="

# Find Prokka location
PROKKA_PATH=$(which prokka 2>/dev/null)

if [ -z "$PROKKA_PATH" ]; then
    echo "ERROR: Prokka not found in PATH"
    echo "Please activate the bactpipe environment first:"
    echo "  micromamba activate bactpipe"
    exit 1
fi

echo "Found Prokka at: $PROKKA_PATH"

# Backup original
if [ ! -f "${PROKKA_PATH}.backup" ]; then
    cp "$PROKKA_PATH" "${PROKKA_PATH}.backup"
    echo "Backup created: ${PROKKA_PATH}.backup"
fi

# Apply the fix - change MINVER from 2.2 to 2.0
sed -i 's/MINVER  => "2.2"/MINVER  => "2.0"/g' "$PROKKA_PATH"

# Verify the fix was applied
if grep -q 'MINVER  => "2.0"' "$PROKKA_PATH"; then
    echo "✓ Fix applied successfully!"
else
    echo "✗ Fix may not have been applied. Check manually."
fi

# Setup Prokka databases
echo "Setting up Prokka databases..."
prokka --setupdb --dbdir ~/prokka_db 2>/dev/null || true

echo "========================================="
echo "Prokka is now fixed and ready to use!"
echo "========================================="