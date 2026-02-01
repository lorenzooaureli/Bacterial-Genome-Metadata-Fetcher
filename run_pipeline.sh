#!/bin/bash

# Bacterial Genome Metadata Fetcher Pipeline
# This script runs the full workflow to retrieve genome accessions and metadata
# from NCBI for bacterial species.

# =============================================================================
# CONFIGURATION
# =============================================================================

# NCBI credentials - set these as environment variables or modify here
# To set as environment variables:
#   export NCBI_EMAIL="your_email@example.com"
#   export NCBI_API_KEY="your_api_key_here"

EMAIL="${NCBI_EMAIL:-}"
API_KEY="${NCBI_API_KEY:-}"

# Number of BioSamples to retrieve per species
BATCH_SIZE="50"

# Check if credentials are set
if [ -z "$EMAIL" ] || [ -z "$API_KEY" ]; then
    echo "Error: NCBI credentials not set."
    echo ""
    echo "Please set your credentials using environment variables:"
    echo "  export NCBI_EMAIL=\"your_email@example.com\""
    echo "  export NCBI_API_KEY=\"your_api_key_here\""
    echo ""
    echo "Or edit this script and set EMAIL and API_KEY directly."
    exit 1
fi

# =============================================================================
# INPUT FILES
# =============================================================================

# List of input .txt files containing species names (one per line)
input_files=(
  environmental.txt
)

# =============================================================================
# PIPELINE EXECUTION
# =============================================================================

for input in "${input_files[@]}"; do
    # Check if input file exists
    if [ ! -f "$input" ]; then
        echo "Warning: Input file '$input' not found. Skipping."
        continue
    fi

    genus=$(basename "$input" .txt)

    echo "=== Processing $genus ==="

    bios_csv="${genus}_bios.csv"
    bios_attr_csv="${genus}_bios_attr.csv"
    bios_gen_csv="${genus}_bios_gen.csv"

    # Step 1: Search for BioSample IDs
    echo "Step 1: Fetching BioSample IDs..."
    python 1_serovars_biosamples_human.py "$input" "$bios_csv" "$EMAIL" "$API_KEY" "$BATCH_SIZE" non-human

    # Step 2: Fetch BioSample metadata
    echo "Step 2: Fetching BioSample metadata..."
    python 2_samp_information_serovars_parallel.py "$bios_csv" "$bios_attr_csv" "$EMAIL" "$API_KEY"

    # Step 3: Fetch genome assembly accessions
    echo "Step 3: Fetching genome accessions..."
    python 3_genomes_serovars_parallel.py "$bios_attr_csv" "$bios_gen_csv" "$EMAIL" "$API_KEY" --filter_gen

    echo "--- Finished $genus ---"
    echo
done

echo "All pipelines completed."
