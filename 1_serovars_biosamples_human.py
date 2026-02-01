# -- Script to search for BioSample IDs of strains corresponding to each organism listed in a text file --

# Example usage: python serovars_biosamples_human.py <input.txt> <output.csv> <email@example.com> <ncbi_apikey> <biosamples_per_species> <non-human>

# Input text file example:

# Vibrio cholerae
# Vibrio harveyi
# vibrio parahaemolyticus

# Example query with non-human flag: "Vibrio harveyi" --> ("Vibrio harveyi"[Organism] OR Vibrio harveyi[All Fields]) NOT (("Homo sapiens"[Host] OR human[Host]) OR ("Homo sapiens"[isolation source] OR human[isolation source]))
# Example query without flag: "Vibrio harveyi" --> ("Vibrio harveyi"[Organism] OR Vibrio harveyi[All Fields])

import pandas as pd
import time
import sys
from Bio import Entrez

# Check if the correct number of arguments is provided
if len(sys.argv) < 5:
    print("Usage: python script.py <input_file> <output_file> <Entrez.email> <Entrez.api_key> <max_results> [non-human]")
    sys.exit(1)

# Read command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
Entrez.email = sys.argv[3]
Entrez.api_key = sys.argv[4]
max_results = int(sys.argv[5])

# Check if 'non-human' is passed as an argument
use_nonhuman_filter = "non-human" in sys.argv

# Configurable retry delay
retry_delay = 5  # Adjust this value as needed (seconds)

def search_biosample(organism_name, max_results=max_results, max_retries=3, use_nonhuman_filter=False):
    """
    Search NCBI BioSample for a given organism name, optionally excluding Homo sapiens samples.
    """
    # Construct the search query
    if use_nonhuman_filter:
        # Exclude samples where host is Homo sapiens or where human appears in host/isolation source
        # Note: Use NOT without AND for correct NCBI boolean logic
        query = f'("{organism_name}"[Organism] OR {organism_name}[All Fields]) NOT (("Homo sapiens"[Host] OR human[Host]) OR ("Homo sapiens"[isolation source] OR human[isolation source]))'
    else:
        query = f'("{organism_name}"[Organism] OR {organism_name}[All Fields])'

    # Print the query for verification
    print(f"🔍 Searching for: {query}")

    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(db="biosample", term=query, retmax=max_results)
            record = Entrez.read(handle)
            handle.close()

            if "IdList" in record and record["IdList"]:
                return [(organism_name, biosample_id) for biosample_id in record["IdList"]]
            else:
                print(f"⚠️ No results found for '{organism_name}'. Query used: {query}")
                return []

        except (IOError, RuntimeError) as e:
            print(f"⚠️ Request failed for '{organism_name}'. Retrying ({attempt + 1}/{max_retries})... Error: {e}")
            time.sleep(retry_delay)

    print(f"❌ Failed to fetch results for '{organism_name}' after {max_retries} attempts.")
    return []

# Read organism names from a text file
try:
    with open(input_file, "r", encoding="utf-8") as file:
        organism_names = [line.strip() for line in file if line.strip()]
except FileNotFoundError:
    print(f"❌ Error: File '{input_file}' not found.")
    sys.exit(1)

# Collect all results
all_results = []
for organism in organism_names:
    all_results.extend(search_biosample(organism, use_nonhuman_filter=use_nonhuman_filter))

# Convert results to a DataFrame
df = pd.DataFrame(all_results, columns=["Organism Name", "BioSample ID"])

# Save the results to a CSV file
df.to_csv(output_file, index=False)

print(f"✅ Results saved to {output_file}")
