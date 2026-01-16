import os
import shutil

def get_base_identifier(filename):
    """
    Extracts 'Gene_EnsemblID' from 'Gene_EnsemblID_fishes.fasta'.
    """
    # 1. Remove extension
    base = os.path.splitext(filename)[0]

    # 2. Split by underscore
    parts = base.split('_')

    # 3. Take the first two parts (Gene and ID)
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    else:
        return parts[0]

def merge_folders(ensembl_dir, homologs_dir, output_dir):
    """
    Matches files between Ensembl and Homolog folders and merges them.
    """
    # 1. Setup Output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # 2. Scan Ensembl Folder (Source of Truth)
    # We iterate over Ensembl files because they have the longer name structure
    ensembl_files = [f for f in os.listdir(ensembl_dir) if f.endswith(('.fasta', '.fa'))]

    if not ensembl_files:
        print(f"No FASTA files found in {ensembl_dir}")
        return

    print(f"Found {len(ensembl_files)} Ensembl files. Starting merge...")
    print("-" * 60)

    count_merged = 0
    count_missing = 0

    for ens_file in ensembl_files:
        # A. Identify the Target Name
        # Input:  ABHD11_ENSG00000106077_fishes.fasta
        # Target: ABHD11_ENSG00000106077.fasta (Expected name in homologs folder)

        identifier = get_base_identifier(ens_file)
        homolog_file = f"{identifier}.fasta"

        # B. Define Full Paths
        path_ens = os.path.join(ensembl_dir, ens_file)
        path_hom = os.path.join(homologs_dir, homolog_file)
        path_out = os.path.join(output_dir, homolog_file) # Final name is the simplified one

        # C. Check if Match Exists
        if os.path.exists(path_hom):
            try:
                with open(path_out, 'w') as outfile:
                    # 1. Write Ensembl Data
                    with open(path_ens, 'r') as f_ens:
                        outfile.write(f_ens.read())
                        # Ensure newline separation
                        if not f_ens.read().endswith('\n'):
                            outfile.write('\n')

                    # 2. Write Homolog Data
                    with open(path_hom, 'r') as f_hom:
                        outfile.write(f_hom.read())
                        # Ensure newline at very end
                        if not f_hom.read().endswith('\n'):
                            outfile.write('\n')

                count_merged += 1

            except Exception as e:
                print(f"  [Error] Failed to merge {identifier}: {e}")

        else:
            # Match not found
            print(f"  [Warning] Missing homolog file for: {identifier}")
            print(f"            (Looked for: {homolog_file})")
            count_missing += 1

    print("-" * 60)
    print(f"Processing Complete.")
    print(f"  Merged:  {count_merged}")
    print(f"  Missing: {count_missing}")
    print(f"  Output:  {output_dir}")

# ==========================================
# CONFIGURATION
# ==========================================
if __name__ == "__main__":

    # 1. Path to your ENSEMBL Orthologs (Input 1)
    # Contains: "ABHD11_ENSG00000106077_fishes.fasta"
    ENSEMBL_FOLDER = "./Downloads"

    # 2. Path to your COMBINED HOMOLOGS (Input 2)
    # Contains: "ABHD11_ENSG00000106077.fasta"
    HOMOLOGS_FOLDER = "./combined_homologs"

    # 3. Path for the FINAL MERGED output
    OUTPUT_FOLDER = "combined_ortho_homologs"

    merge_folders(ENSEMBL_FOLDER, HOMOLOGS_FOLDER, OUTPUT_FOLDER)
