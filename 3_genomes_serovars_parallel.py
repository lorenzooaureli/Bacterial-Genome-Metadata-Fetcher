# -- Script to fetch genome assembly information for BioSamples in parallel --
#
# Usage options:
# 1. Positional arguments:
#    python genomes_serovars_parallel.py input.csv output.csv email@example.com api_key [num_batches (default=5)]
#
# 2. Named arguments:
#    python genomes_serovars_parallel.py --input_file input.csv --output_file output.csv --email email@example.com --api_key api_key [--num_batches 10]
#
# Additional options:
#    --retries 5           Number of retries for failed requests (default: 5)
#    --delay 10            Delay between retries in seconds (default: 10)
#    --rate_limit 3        Maximum requests per second to NCBI (default: 3)
#    --filter_gen          Filter out rows with no assembly information
#                          (keeps only rows with values in either 'NCBI Assembly Accession' or 'RefSeq Assembly (GCF_)')
#
# Example:
#    python genomes_serovars_parallel.py K_aerogenes_biosamp.csv K_aerogenes_genomes.csv laureli@lbl.gov 12e3b20bbfc7b01f7892f2421fbfd0370d09 10 --filter_gen

import pandas as pd
from Bio import Entrez
import time
import concurrent.futures
import argparse
import sys
from tqdm import tqdm
import numpy as np
import random
import threading

# Create a rate limiter class to control API requests
class RateLimiter:
    def __init__(self, max_per_second):
        self.lock = threading.Lock()
        self.max_per_second = max_per_second
        self.tokens = max_per_second
        self.updated_at = time.monotonic()

    def acquire(self):
        """Acquire a token to make an API request"""
        with self.lock:
            now = time.monotonic()
            # Add tokens based on time elapsed, up to max_per_second
            elapsed = now - self.updated_at
            new_tokens = elapsed * self.max_per_second
            if new_tokens > 0:
                self.tokens = min(self.tokens + new_tokens, self.max_per_second)
                self.updated_at = now

            # If we have tokens, consume one and proceed
            if self.tokens >= 1:
                self.tokens -= 1
                return 0

            # If no tokens, calculate sleep time and wait
            sleep_time = (1 - self.tokens) / self.max_per_second
            return sleep_time

    def wait(self):
        """Wait until a token is available"""
        sleep_time = self.acquire()
        if sleep_time > 0:
            time.sleep(sleep_time)
            # Add a small random delay to prevent threads from synchronizing
            time.sleep(random.uniform(0, 0.1))

def parse_arguments():
    """Parse command line arguments, supporting both positional and named arguments."""
    parser = argparse.ArgumentParser(description='Fetch assembly details for BioSamples in parallel.')

    # Create argument groups to handle both positional and named arguments
    positional_group = parser.add_argument_group('Positional Arguments')
    named_group = parser.add_argument_group('Named Arguments')

    # Add positional arguments
    positional_group.add_argument('input_file', nargs='?', type=str, help='Input CSV file with BioSample IDs')
    positional_group.add_argument('output_file', nargs='?', type=str, help='Output CSV file for results')
    positional_group.add_argument('email', nargs='?', type=str, help='Email for NCBI API')
    positional_group.add_argument('api_key', nargs='?', type=str, help='NCBI API key')
    positional_group.add_argument('num_batches', nargs='?', type=int, help='Number of parallel batches to process')

    # Add named arguments (with -- prefix)
    named_group.add_argument('--input_file', dest='named_input_file', type=str, help='Input CSV file with BioSample IDs')
    named_group.add_argument('--output_file', dest='named_output_file', type=str, help='Output CSV file for results')
    named_group.add_argument('--email', dest='named_email', type=str, help='Email for NCBI API')
    named_group.add_argument('--api_key', dest='named_api_key', type=str, help='NCBI API key')
    named_group.add_argument('--num_batches', dest='named_num_batches', type=int, help='Number of parallel batches to process')

    # Add other optional arguments
    parser.add_argument('--retries', type=int, default=5, help='Number of retries for failed requests')
    parser.add_argument('--delay', type=float, default=10, help='Delay between retries in seconds')
    parser.add_argument('--rate_limit', type=float, default=3, help='Maximum requests per second to NCBI')
    parser.add_argument('--filter_gen', action='store_true', help='Filter out rows with no assembly information')

    args = parser.parse_args()

    # Resolve between positional and named arguments
    resolved_args = argparse.Namespace()

    # For each parameter, use the named version if provided, otherwise use the positional version
    resolved_args.input_file = args.named_input_file if args.named_input_file else args.input_file
    resolved_args.output_file = args.named_output_file if args.named_output_file else args.output_file
    resolved_args.email = args.named_email if args.named_email else args.email
    resolved_args.api_key = args.named_api_key if args.named_api_key else args.api_key

    # For num_batches, use named version if provided, otherwise positional, otherwise default to 10
    if args.named_num_batches is not None:
        resolved_args.num_batches = args.named_num_batches
    elif args.num_batches is not None:
        resolved_args.num_batches = args.num_batches
    else:
        resolved_args.num_batches = 5

    # Copy over other arguments
    resolved_args.retries = args.retries
    resolved_args.delay = args.delay
    resolved_args.rate_limit = args.rate_limit
    resolved_args.filter_gen = args.filter_gen

    # Validate required arguments
    if not resolved_args.input_file or not resolved_args.output_file or not resolved_args.email or not resolved_args.api_key:
        parser.error("Missing required arguments. Please provide input_file, output_file, email, and api_key.")

    return resolved_args

def get_assembly_details(biosample_id, email, api_key, retries=5, delay=10, rate_limiter=None):
    """
    Get assembly details for a BioSample ID.

    Args:
        biosample_id: The BioSample ID to query
        email: Email for NCBI API
        api_key: NCBI API key
        retries: Number of retries for failed requests
        delay: Delay between retries in seconds
        rate_limiter: Rate limiter instance to control API request rate

    Returns:
        Tuple of (biosample_id, ncbi_assembly, refseq)
    """
    # Set NCBI credentials for this thread
    Entrez.email = email
    Entrez.api_key = api_key

    for attempt in range(retries):
        try:
            # Wait for rate limiter if provided
            if rate_limiter:
                rate_limiter.wait()

            # Search the Assembly database using the BioSample ID
            search_handle = Entrez.esearch(db="assembly", term=biosample_id, retmode="xml")
            search_record = Entrez.read(search_handle)
            search_handle.close()

            if not search_record["IdList"]:
                return biosample_id, None, None  # No results

            assembly_uid = search_record["IdList"][0]  # Get first assembly UID

            # Wait for rate limiter again before second request
            if rate_limiter:
                rate_limiter.wait()

            # Fetch Assembly details
            fetch_handle = Entrez.esummary(db="assembly", id=assembly_uid, retmode="xml")
            summary_record = Entrez.read(fetch_handle)
            fetch_handle.close()

            # Extract assembly accessions
            doc = summary_record['DocumentSummarySet']['DocumentSummary'][0]
            ncbi_assembly = doc.get("AssemblyAccession", None)  # GCA_xxxxxx or PDT_xxxx
            refseq = doc.get("RefSeq", None)  # GCF_xxxxxx (if available)

            return biosample_id, ncbi_assembly, refseq

        except Exception as e:
            if attempt < retries - 1:  # Don't sleep on the last attempt
                print(f"Attempt {attempt + 1} failed for {biosample_id}: {e}")
                time.sleep(delay + random.uniform(0, 1))  # Add jitter to delay
            else:
                print(f"All attempts failed for {biosample_id}: {e}")

    # If all attempts fail, return None values
    return biosample_id, None, None

def process_batch(batch, email, api_key, retries, delay, rate_limit):
    """Process a batch of BioSample IDs."""
    # Create a rate limiter for this batch
    rate_limiter = RateLimiter(rate_limit)

    results = []
    for biosample in batch:
        results.append(get_assembly_details(biosample, email, api_key, retries, delay, rate_limiter))
    return results

def main():
    # Parse command line arguments
    args = parse_arguments()

    # Load dataset
    print(f"Loading data from {args.input_file}")
    try:
        df = pd.read_csv(args.input_file)
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

    # Ensure df has the "BioSample" column
    if "BioSample" not in df.columns:
        print("Error: The input file does not contain a 'BioSample' column.")
        sys.exit(1)

    # Extract unique BioSample IDs
    biosample_ids = df["BioSample"].dropna().unique().tolist()
    print(f"Found {len(biosample_ids)} unique BioSample IDs")

    # Split BioSample IDs into batches
    num_batches = min(args.num_batches, len(biosample_ids))
    batches = np.array_split(biosample_ids, num_batches)
    print(f"Processing in {num_batches} parallel batches")

    # Process batches in parallel
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_batches) as executor:
        # Submit all batches to the executor
        future_to_batch = {
            executor.submit(process_batch, list(batch), args.email, args.api_key, args.retries, args.delay, args.rate_limit): i
            for i, batch in enumerate(batches)
        }

        # Process results as they complete
        for future in tqdm(concurrent.futures.as_completed(future_to_batch),
                          total=len(future_to_batch),
                          desc="Processing batches"):
            batch_index = future_to_batch[future]
            try:
                batch_results = future.result()
                results.extend(batch_results)
                print(f"Batch {batch_index+1}/{num_batches} completed")
            except Exception as e:
                print(f"Batch {batch_index+1}/{num_batches} failed: {str(e)}")

    # Convert results to DataFrame
    print("Creating results dataframe")
    df_results = pd.DataFrame(results, columns=["BioSample", "NCBI Assembly Accession", "RefSeq Assembly (GCF_)"])

    # Merge with original dataset
    df_final = df.merge(df_results, on="BioSample", how="left")

    # Filter out rows with no assembly information if requested
    if args.filter_gen:
        print("Filtering out rows with no assembly information...")

        # Count rows before filtering
        rows_before = len(df_final)

        # Create a function to check if a value is not empty
        def has_value(value):
            return not pd.isna(value) and value is not None and str(value).strip() != ""

        # Filter the dataframe
        # Keep rows where at least one of the assembly columns has a value
        filtered_df = df_final[
            df_final["NCBI Assembly Accession"].apply(has_value) |
            df_final["RefSeq Assembly (GCF_)"].apply(has_value)
        ]

        # Count rows after filtering
        rows_after = len(filtered_df)
        rows_removed = rows_before - rows_after

        print(f"Removed {rows_removed} rows with no assembly information.")
        print(f"Kept {rows_after} rows with at least one assembly accession.")

        # Use the filtered dataframe
        df_final = filtered_df

    # Save results
    print(f"Saving results to {args.output_file}")
    df_final.to_csv(args.output_file, index=False)
    print("Done!")

if __name__ == "__main__":
    main()
