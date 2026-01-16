#!/usr/bin/env python3
import requests
import sys
import time
import pandas as pd
import os

# --- CONFIGURATION ---
SERVER = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}
OUTPUT_DIR = "Downloads"
UNIQUE_LIST_FILENAME = "unique_gene_list.txt"

# --- TAXONOMY LEVEL FILTER ---
# These are the ancestral nodes that contain fishes but exclude Tetrapods (Mammals/Birds).
FISH_TAXONOMY_LEVELS = {
    "Vertebrata",
    "Gnathostomata",
    "Euteleostomi",
    "Sarcopterygii",
    "Actinopterygii"
}

def fetch_url(endpoint, params=None):
    """
    Fetches data with robust retry logic for rate limits (HTTP 429).
    """
    url = SERVER + endpoint
    retries = 5
    for attempt in range(retries):
        try:
            r = requests.get(url, headers=HEADERS, params=params, timeout=60)
            if r.status_code == 429:
                wait = float(r.headers.get("Retry-After", 2))
                time.sleep(wait)
                continue
            if r.ok: return r.json()
            else:
                if r.status_code == 404: return None
                return None
        except requests.exceptions.RequestException:
            time.sleep(2)
    return None

def get_gene_symbol(gene_id):
    endpoint = f"/lookup/id/{gene_id}"
    data = fetch_url(endpoint)
    if data and 'display_name' in data:
        return data['display_name']
    return gene_id

def get_orthologs(gene_id):
    """
    Fetches orthologs and filters STRICTLY by the taxonomy_level field.
    """
    endpoint = f"/homology/id/human/{gene_id}"
    params = {"type": "orthologues", "format": "condensed"}

    data = fetch_url(endpoint, params)
    results = []

    if not data or 'data' not in data or not isinstance(data['data'], list) or len(data['data']) == 0:
        return results

    homologies = data['data'][0].get('homologies', [])

    for h in homologies:
        tax_level = h.get('taxonomy_level')

        # Filter: Only keep if the split occurred at a "Fish" node
        if tax_level in FISH_TAXONOMY_LEVELS:
            results.append({
                'source_gene': gene_id,
                'species': h.get('species'),
                'target_gene_id': h.get('id'),
                'target_protein_id': h.get('protein_id'),
                'homology_type': h.get('type'),
                'taxonomy_level': tax_level
            })

    return results

def get_transcript_and_cds(protein_id):
    if not protein_id: return None, None

    # 1. Get Parent Transcript
    endpoint_trans = f"/overlap/translation/{protein_id}"
    data_trans = fetch_url(endpoint_trans)

    transcript_id = None
    if data_trans and len(data_trans) > 0:
        transcript_id = data_trans[0].get('Parent')

    if not transcript_id: return None, None

    # 2. Get CDS Sequence
    endpoint_seq = f"/sequence/id/{transcript_id}"
    data_seq = fetch_url(endpoint_seq, params={"type": "cds"})

    sequence = None
    if data_seq and 'seq' in data_seq:
        sequence = data_seq['seq']
    return transcript_id, sequence

def create_unique_list(input_file):
    """
    Reads the input file, removes duplicates while preserving order,
    and writes the result to a new file.
    """
    unique_genes = []
    seen = set()

    # Read and Deduplicate
    with open(input_file, 'r') as f:
        for line in f:
            g = line.strip()
            if g and g not in seen:
                seen.add(g)
                unique_genes.append(g)

    # Write to new file
    with open(UNIQUE_LIST_FILENAME, 'w') as f:
        for g in unique_genes:
            f.write(g + "\n")

    print(f"--- PRE-PROCESSING ---")
    print(f"Input file: {input_file}")
    print(f"Duplicates removed. Unique list saved to: {UNIQUE_LIST_FILENAME}")
    print(f"Total Unique Genes to Process: {len(unique_genes)}")
    print(f"----------------------\n")

    return unique_genes

def main():
    # --- 1. SETUP ---
    if len(sys.argv) < 2:
        print("Usage: python3 fetch_fishes_log_unique.py <gene_ids.txt>")
        sys.exit(1)

    input_arg = sys.argv[1]

    # Only proceed if it's a file
    if os.path.exists(input_arg):
        unique_genes = create_unique_list(input_arg)
    else:
        # Fallback for single ID testing
        unique_genes = [input_arg]
        print(f"Processing single ID input: {input_arg}")

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    print(f"Starting Download...")
    print(f"Filtering for Taxonomy Levels: {FISH_TAXONOMY_LEVELS}")

    # --- 2. MAIN LOOP ---
    for i, gene in enumerate(unique_genes):
        print(f"[{i+1}/{len(unique_genes)}] Processing {gene}...")

        gene_name = get_gene_symbol(gene)
        print(f"    -> Identified as: {gene_name}")

        # Filenames
        csv_filename = os.path.join(OUTPUT_DIR, f"{gene_name}_{gene}_fishes.csv")
        fasta_filename = os.path.join(OUTPUT_DIR, f"{gene_name}_{gene}_fishes.fasta")

        # Fetch Orthologs
        orthologs = get_orthologs(gene)

        if not orthologs:
            print(f"    -> No matching fish orthologs found.")
            continue

        print(f"    -> Found {len(orthologs)} fish matches.")

        gene_metadata = []
        gene_fasta_lines = []
        count = 0

        # Fetch Sequences
        for ortho in orthologs:
            pid = ortho['target_protein_id']
            if not pid: continue

            tid, seq = get_transcript_and_cds(pid)

            if tid and seq:
                ortho['transcript_id'] = tid
                ortho['source_gene_name'] = gene_name
                gene_metadata.append(ortho)

                # Header format: >Transcript | Protein | Species | GeneName | GeneID | TaxLevel
                header = f">{tid} | {pid} | {ortho['species']} | {gene_name} | {gene} | {ortho['taxonomy_level']}"
                gene_fasta_lines.append(f"{header}\n{seq}")
                count += 1

        # Save Files
        if gene_metadata:
            df = pd.DataFrame(gene_metadata)
            # Reorder columns
            cols = ['source_gene_name', 'source_gene', 'species', 'taxonomy_level'] + [c for c in df.columns if c not in ['source_gene_name', 'source_gene', 'species', 'taxonomy_level']]
            df = df[cols]
            df.to_csv(csv_filename, index=False)

            with open(fasta_filename, "w") as f:
                f.write("\n".join(gene_fasta_lines))
            print(f"    -> Saved: {csv_filename}")
            print(f"    -> Saved: {fasta_filename} ({count} seqs)")
        else:
            print("    -> No valid CDS sequences retrieved.")

    print("\nAll Done.")

if __name__ == "__main__":
    main()
