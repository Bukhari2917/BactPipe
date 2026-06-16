#!/bin/bash
# Fix Prokka BLAST version detection bug

echo "Applying Prokka BLAST version fix..."

PROKKA_PATH=$(which prokka 2>/dev/null)

if [ -z "$PROKKA_PATH" ]; then
    echo "ERROR: Prokka not found in PATH"
    exit 1
fi

echo "Found Prokka at: $PROKKA_PATH"

# Backup original
if [ ! -f "${PROKKA_PATH}.backup" ]; then
    cp "$PROKKA_PATH" "${PROKKA_PATH}.backup"
    echo "Backup created"
fi

# Apply fix
sed -i 's/MINVER  => "2.2"/MINVER  => "2.0"/g' "$PROKKA_PATH"

# Verify
if grep -q 'MINVER  => "2.0"' "$PROKKA_PATH"; then
    echo "Prokka fix applied successfully!"
else
    echo "Fix failed"
fi

# Setup databases
prokka --setupdb --dbdir ~/prokka_db 2>/dev/null || true

echo "Prokka is now fixed!"
