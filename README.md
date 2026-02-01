# Bacterial Genome Metadata Fetcher

A pipeline to retrieve genome accessions and metadata from NCBI for bacterial species. Starting from a list of bacterial species names, this workflow fetches BioSample information, detailed metadata (isolation source, host, serovar, etc.), and associated genome assembly accessions.

## Overview

This pipeline consists of three Python scripts that run sequentially:

1. **Step 1**: Search NCBI BioSample database for each species and retrieve BioSample IDs
2. **Step 2**: Fetch detailed metadata for each BioSample (isolation source, host, strain, serovar, etc.)
3. **Step 3**: Retrieve associated genome assembly accessions (GCA/GCF) for each BioSample

## Requirements

### Dependencies

```bash
pip install biopython pandas numpy tqdm
```

### NCBI Credentials

You need an NCBI account and API key to use this pipeline:

1. Create an NCBI account at https://www.ncbi.nlm.nih.gov/account/
2. Go to your account settings and generate an API key
3. Set your credentials as environment variables:

```bash
export NCBI_EMAIL="your_email@example.com"
export NCBI_API_KEY="your_api_key_here"
```

Or create a `.env` file (do not commit this to version control):

```
NCBI_EMAIL=your_email@example.com
NCBI_API_KEY=your_api_key_here
```

## Input Format

Create a text file with one bacterial species name per line:

```
Deinococcus radiodurans
Bacillus subtilis
Escherichia coli
```

## Usage

### Running the Full Pipeline

```bash
# Set your credentials
export NCBI_EMAIL="your_email@example.com"
export NCBI_API_KEY="your_api_key_here"

# Run the pipeline
./run_pipeline.sh
```

### Running Scripts Individually

Each script can be run independently with either positional or named arguments.

---

## Script 1: `1_serovars_biosamples_human.py`

Searches the NCBI BioSample database for each species name and retrieves BioSample IDs.

### Usage

```bash
python 1_serovars_biosamples_human.py <input.txt> <output.csv> <email> <api_key> <max_results> [non-human]
```

### Arguments

| Argument | Description |
|----------|-------------|
| `input.txt` | Text file with species names (one per line) |
| `output.csv` | Output CSV file for BioSample IDs |
| `email` | Your NCBI registered email |
| `api_key` | Your NCBI API key |
| `max_results` | Maximum number of BioSamples to retrieve per species |
| `non-human` | Optional flag to exclude samples from human hosts |

### Flags

- **`non-human`**: When specified, excludes BioSamples where the host is *Homo sapiens* or human. This is useful for environmental or non-clinical studies.

### Example

```bash
# Include all samples
python 1_serovars_biosamples_human.py species_list.txt biosamples.csv $NCBI_EMAIL $NCBI_API_KEY 50

# Exclude human-derived samples
python 1_serovars_biosamples_human.py species_list.txt biosamples.csv $NCBI_EMAIL $NCBI_API_KEY 50 non-human
```

### Output

CSV file with columns:
- `Organism Name`: The queried species name
- `BioSample ID`: NCBI BioSample identifier

---

## Script 2: `2_samp_information_serovars_parallel.py`

Fetches detailed metadata for each BioSample ID using parallel processing.

### Usage

**Positional arguments:**
```bash
python 2_samp_information_serovars_parallel.py <input.csv> <output.csv> <email> <api_key> [num_batches]
```

**Named arguments:**
```bash
python 2_samp_information_serovars_parallel.py --input <input.csv> --output <output.csv> --email <email> --api_key <api_key> [options]
```

### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `input.csv` | Input CSV with BioSample IDs (from Step 1) | Required |
| `output.csv` | Output CSV for metadata | Required |
| `email` | Your NCBI registered email | Required |
| `api_key` | Your NCBI API key | Required |
| `num_batches` / `--batches` | Number of parallel batches | 5 |

### Optional Flags

| Flag | Description | Default |
|------|-------------|---------|
| `--retries` | Number of retries for failed API requests | 5 |
| `--delay` | Delay in seconds between retries | 10 |
| `--rate_limit` | Maximum API requests per second | 3 |
| `--filter_empty` | Remove rows with no meaningful values in Source Type, Host Disease, or Isolation Source | Off |

### Example

```bash
# Basic usage
python 2_samp_information_serovars_parallel.py biosamples.csv metadata.csv $NCBI_EMAIL $NCBI_API_KEY

# With filtering and more batches
python 2_samp_information_serovars_parallel.py biosamples.csv metadata.csv $NCBI_EMAIL $NCBI_API_KEY 10 --filter_empty

# Using named arguments
python 2_samp_information_serovars_parallel.py \
    --input biosamples.csv \
    --output metadata.csv \
    --email $NCBI_EMAIL \
    --api_key $NCBI_API_KEY \
    --batches 10 \
    --rate_limit 5 \
    --filter_empty
```

### Output

CSV file with columns:
- `BioSample`: BioSample accession
- `Species`: Organism taxonomy name
- `Strain`: Strain identifier
- `Isolation Source`: Where the sample was isolated from
- `Source Type`: Type of source (e.g., environmental, clinical)
- `Host`: Host organism
- `Host Disease`: Disease associated with the host
- `Serovar`: Serovar designation
- `Serogroup`: Serogroup designation
- `Serotype`: Serotype designation
- `Serovar or Serotype`: Combined serovar/serotype field

---

## Script 3: `3_genomes_serovars_parallel.py`

Retrieves genome assembly accessions for each BioSample using parallel processing.

### Usage

**Positional arguments:**
```bash
python 3_genomes_serovars_parallel.py <input.csv> <output.csv> <email> <api_key> [num_batches]
```

**Named arguments:**
```bash
python 3_genomes_serovars_parallel.py --input_file <input.csv> --output_file <output.csv> --email <email> --api_key <api_key> [options]
```

### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `input.csv` | Input CSV with BioSample data (from Step 2) | Required |
| `output.csv` | Output CSV with genome accessions | Required |
| `email` | Your NCBI registered email | Required |
| `api_key` | Your NCBI API key | Required |
| `num_batches` / `--num_batches` | Number of parallel batches | 5 |

### Optional Flags

| Flag | Description | Default |
|------|-------------|---------|
| `--retries` | Number of retries for failed API requests | 5 |
| `--delay` | Delay in seconds between retries | 10 |
| `--rate_limit` | Maximum API requests per second | 3 |
| `--filter_gen` | Remove rows with no genome assembly information | Off |

### Example

```bash
# Basic usage
python 3_genomes_serovars_parallel.py metadata.csv genomes.csv $NCBI_EMAIL $NCBI_API_KEY

# With genome filtering
python 3_genomes_serovars_parallel.py metadata.csv genomes.csv $NCBI_EMAIL $NCBI_API_KEY 10 --filter_gen

# Using named arguments
python 3_genomes_serovars_parallel.py \
    --input_file metadata.csv \
    --output_file genomes.csv \
    --email $NCBI_EMAIL \
    --api_key $NCBI_API_KEY \
    --num_batches 10 \
    --rate_limit 5 \
    --filter_gen
```

### Output

The input CSV is augmented with additional columns:
- `NCBI Assembly Accession`: GenBank assembly accession (GCA_xxxxxx)
- `RefSeq Assembly (GCF_)`: RefSeq assembly accession (GCF_xxxxxx), if available

---

## Complete Workflow Example

```bash
# 1. Set credentials
export NCBI_EMAIL="your_email@example.com"
export NCBI_API_KEY="your_api_key_here"

# 2. Create input file
echo -e "Bacillus subtilis\nEscherichia coli\nPseudomonas aeruginosa" > species.txt

# 3. Run the pipeline
python 1_serovars_biosamples_human.py species.txt biosamples.csv $NCBI_EMAIL $NCBI_API_KEY 50 non-human
python 2_samp_information_serovars_parallel.py biosamples.csv metadata.csv $NCBI_EMAIL $NCBI_API_KEY
python 3_genomes_serovars_parallel.py metadata.csv final_output.csv $NCBI_EMAIL $NCBI_API_KEY --filter_gen

# 4. Results are in final_output.csv
```

## Notes

- The NCBI API has rate limits. The scripts include built-in rate limiting (default: 3 requests/second) and retry logic.
- Using an API key increases your rate limit from 3 to 10 requests per second.
- For large datasets, consider increasing the number of parallel batches to speed up processing.
- The `--filter_*` flags are useful for cleaning the output to include only samples with relevant data.

## License

MIT License
