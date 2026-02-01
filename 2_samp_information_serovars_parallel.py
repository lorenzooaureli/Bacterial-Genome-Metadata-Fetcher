# -- Script to gather genomes and biosample information from biosample accession list in parallel --
# Usage options:
# 1. Positional arguments:
#    python samp_information_serovars_parallel.py input.csv output.csv email@example.com api_key [num_batches (default=5)]
#
# 2. Named arguments:
#    python samp_information_serovars_parallel.py --input input.csv --output output.csv --email email@example.com --api_key api_key [--batches 10]
#
# Additional options:
#    --retries 5           Number of retries for failed requests (default: 5)
#    --delay 10            Delay between retries in seconds (default: 10)
#    --rate_limit 3        Maximum requests per second to NCBI (default: 3)
#    --filter_empty        Filter out rows with no meaningful values in key columns
#                          (Source Type, Host Disease, Isolation Source)
#
# Example:
#    python samp_information_serovars_parallel.py K_aerogenes.csv K_aerogenes_biosamp.csv laureli@lbl.gov 12e3b20bbfc7b01f7892f2421fbfd0370d09 10 --filter_empty
#
# input.csv structure:
# Organism Name	BioSample ID
# Salmonella enterica subsp. arizonae serovar 1,13,23:g,z51:-	12647575
# Salmonella enterica subsp. arizonae serovar 1,13,23:g,z51:-	12642787
# ...

import pandas as pd
import sys
import argparse
from Bio import Entrez
import time
import xml.etree.ElementTree as ET
import concurrent.futures
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
    parser = argparse.ArgumentParser(description='Fetch BioSample details in parallel.')

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
    named_group.add_argument('--input', dest='named_input_file', type=str, help='Input CSV file with BioSample IDs')
    named_group.add_argument('--output', dest='named_output_file', type=str, help='Output CSV file for results')
    named_group.add_argument('--email', dest='named_email', type=str, help='Email for NCBI API')
    named_group.add_argument('--api_key', dest='named_api_key', type=str, help='NCBI API key')
    named_group.add_argument('--batches', dest='named_num_batches', type=int, help='Number of parallel batches to process')

    # Add other optional arguments
    parser.add_argument('--retries', type=int, default=5, help='Number of retries for failed requests')
    parser.add_argument('--delay', type=float, default=10, help='Delay between retries in seconds')
    parser.add_argument('--rate_limit', type=float, default=3, help='Maximum requests per second to NCBI')
    parser.add_argument('--filter_empty', action='store_true', help='Filter out rows with no meaningful values in key columns')

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
    resolved_args.filter_empty = args.filter_empty

    # Validate required arguments
    if not resolved_args.input_file or not resolved_args.output_file or not resolved_args.email or not resolved_args.api_key:
        parser.error("Missing required arguments. Please provide input_file, output_file, email, and api_key.")

    return resolved_args

def fetch_biosample_details_and_assemblies(biosample_accession, email, api_key, retries=5, delay=10, rate_limiter=None):
    """
    Fetch details for a BioSample accession.

    Args:
        biosample_accession: The BioSample ID to query
        email: Email for NCBI API
        api_key: NCBI API key
        retries: Number of retries for failed requests
        delay: Delay between retries in seconds
        rate_limiter: Rate limiter instance to control API request rate

    Returns:
        Dictionary containing BioSample details
    """
    # Set NCBI credentials for this thread
    Entrez.email = email
    Entrez.api_key = api_key

    biosample_details = {
        "BioSample": biosample_accession,
        "Species": "N/A",
        "Strain": "N/A",
        "Isolation Source": "N/A",
        "Source Type": "N/A",
        "Host": "N/A",
        "Host Disease": "N/A",
        "Serovar": "N/A",
        "Serogroup": "N/A",
        "Serotype": "N/A",
        "Serovar or Serotype": "N/A",
        "Assemblies": [],
        "Assembly Levels": []
    }

    for attempt in range(retries):
        try:
            # Wait for rate limiter if provided
            if rate_limiter:
                rate_limiter.wait()

            # Fetch BioSample details
            handle = Entrez.efetch(db="biosample", id=biosample_accession, retmode="xml")
            xml_data = handle.read()
            handle.close()

            # Parse XML using ElementTree
            root = ET.fromstring(xml_data)

            # Extract species name
            organism = root.find(".//Organism")
            if organism is not None:
                biosample_details["Species"] = organism.attrib.get("taxonomy_name", "N/A")

            # Extract attributes
            attributes = root.findall(".//Attribute")
            for attribute in attributes:
                name = attribute.get("harmonized_name")
                value = attribute.text
                if attribute.attrib.get("attribute_name", "").lower() == "strain":
                    biosample_details["Strain"] = value
                if name == "isolation_source":
                    biosample_details["Isolation Source"] = value
                elif name == "source_type":
                    biosample_details["Source Type"] = value
                elif name == "host":
                    biosample_details["Host"] = value
                elif name == "host_disease":
                    biosample_details["Host Disease"] = value
                elif name == "serovar":
                    biosample_details["Serovar"] = value
                elif name == "serogroup":
                    biosample_details["Serogroup"] = value
                elif name == "serotype":
                    biosample_details["Serotype"] = value
                elif name in ["serovar_or_serotype", "serovar or serotype"]:
                    biosample_details["Serovar or Serotype"] = value

            return biosample_details

        except Exception as e:
            if attempt < retries - 1:  # Don't sleep on the last attempt
                print(f"Attempt {attempt + 1} failed for {biosample_accession}: {e}")
                time.sleep(delay)
            else:
                print(f"All attempts failed for {biosample_accession}: {e}")

    return biosample_details  # Return default values if all attempts fail

def process_batch(batch, email, api_key, retries, delay, rate_limit):
    """Process a batch of BioSample IDs."""
    # Create a rate limiter for this batch
    rate_limiter = RateLimiter(rate_limit)

    results = []
    for biosample_id in batch:
        results.append(fetch_biosample_details_and_assemblies(biosample_id, email, api_key, retries, delay, rate_limiter))
    return results

def main():
    # Parse command line arguments
    args = parse_arguments()

    # Read dataset
    print(f"Loading data from {args.input_file}")
    try:
        df = pd.read_csv(args.input_file)
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

    # Ensure df has the "BioSample ID" column
    if "BioSample ID" not in df.columns:
        print("Error: The input file does not contain a 'BioSample ID' column.")
        sys.exit(1)

    # Extract BioSample IDs and add "SAMN" prefix if missing
    biosample_ids = df["BioSample ID"].astype(str)
    biosample_ids = [f"SAMN{id}" if not id.startswith("SAMN") else id for id in biosample_ids]
    print(f"Found {len(biosample_ids)} BioSample IDs")

    # Split BioSample IDs into batches
    num_batches = min(args.num_batches, len(biosample_ids))
    batches = np.array_split(biosample_ids, num_batches)
    print(f"Processing in {num_batches} parallel batches")

    # Process batches in parallel
    all_biosamples = []
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
                all_biosamples.extend(batch_results)
                print(f"Batch {batch_index+1}/{num_batches} completed")
            except Exception as e:
                print(f"Batch {batch_index+1}/{num_batches} failed: {str(e)}")

    # Convert results to DataFrame
    print("Creating results dataframe")
    final_df = pd.DataFrame(all_biosamples)

    # Convert lists to comma-separated strings
    final_df["Assemblies"] = final_df["Assemblies"].apply(lambda x: ", ".join(x) if isinstance(x, list) else x)
    final_df["Assembly Levels"] = final_df["Assembly Levels"].apply(lambda x: ", ".join(x) if isinstance(x, list) else x)

    # Filter out rows with no meaningful values in key columns if requested
    if args.filter_empty:
        print("Filtering out rows with no meaningful values in key columns...")

        # Define the columns to check
        key_columns = ['Source Type', 'Host Disease', 'Isolation Source']

        # Define values to consider as empty/missing
        empty_values = ['missing', 'not collected', 'not provided', 'n/a', '']

        # Count rows before filtering
        rows_before = len(final_df)

        # Create a function to check if a value is meaningful
        def is_meaningful(value):
            if pd.isna(value):
                return False
            if str(value).lower() in empty_values:
                return False
            return True

        # Filter the dataframe
        # Keep rows where at least one of the key columns has a meaningful value
        filtered_df = final_df[final_df[key_columns].applymap(is_meaningful).any(axis=1)]

        # Count rows after filtering
        rows_after = len(filtered_df)
        rows_removed = rows_before - rows_after

        print(f"Removed {rows_removed} rows with no meaningful values in key columns.")
        print(f"Kept {rows_after} rows with at least one meaningful value.")

        # Use the filtered dataframe
        final_df = filtered_df

    # Save results to output file
    if not final_df.empty:
        print(f"Saving results to {args.output_file}")
        final_df.to_csv(args.output_file, index=False)
        print(f"Results saved to {args.output_file}")
    else:
        print("Warning: No data to write.")

if __name__ == "__main__":
    main()
