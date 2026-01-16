import os
import subprocess
import shutil
import sys

# ==========================================
# 1. HELPER FUNCTIONS
# ==========================================

def get_alignment_length(fasta_path):
    """
    Reads the first sequence in the FASTA file to determine alignment length.
    CORRECTION: Removes spaces (formatting) but counts dashes (gaps) and letters.
    """
    seq_len = 0
    with open(fasta_path, 'r') as f:
        header_seen = False
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header_seen:
                    # We reached the next sequence; stop reading.
                    break
                header_seen = True
            else:
                if header_seen:
                    # Remove spaces used for visual formatting (e.g. "gatc gatc")
                    # But KEEP dashes ("-") because they are part of the alignment length
                    clean_line = line.replace(" ", "")
                    seq_len += len(clean_line)
    return seq_len

def create_partition_file(file_path, length):
    """
    Creates a NEXUS/RAxML style partition file for 3 codon positions.
    Format:
    DNA, pos1 = 1-Length\3
    DNA, pos2 = 2-Length\3
    DNA, pos3 = 3-Length\3
    """
    content = (
        f"DNA, pos1 = 1-{length}\\3\n"
        f"DNA, pos2 = 2-{length}\\3\n"
        f"DNA, pos3 = 3-{length}\\3\n"
    )

    with open(file_path, 'w') as f:
        f.write(content)
    return True

def get_gene_id(filename):
    """
    Extracts 'Gene_EnsemblID' from filename.
    Example: 'ABHD11_ENSG00000106077.fasta' -> 'ABHD11_ENSG00000106077'
    """
    base = os.path.splitext(filename)[0]
    parts = base.split('_')
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return parts[0]

# ==========================================
# 2. MAIN PIPELINE
# ==========================================

def run_phylogeny_pipeline(input_folder, output_root):

    # 1. Setup Main Directory
    if not os.path.exists(output_root):
        os.makedirs(output_root)
        print(f"Created main directory: {output_root}")

    # 2. Find Input Files
    fasta_files = [f for f in os.listdir(input_folder) if f.endswith(('.fasta', '.fa'))]

    if not fasta_files:
        print(f"No fasta files found in {input_folder}")
        return

    print(f"Found {len(fasta_files)} alignments. Starting IQ-TREE pipeline...")
    print("=" * 60)

    for filename in fasta_files:
        original_path = os.path.join(input_folder, filename)

        # A. Create Gene-Specific Folder
        gene_id = get_gene_id(filename)
        gene_folder = os.path.join(output_root, gene_id)

        if not os.path.exists(gene_folder):
            os.makedirs(gene_folder)

        # B. Copy Alignment to Gene Folder
        dest_fasta_path = os.path.join(gene_folder, filename)
        shutil.copy2(original_path, dest_fasta_path)

        # C. Get Length and Write Partition File
        aln_len = get_alignment_length(dest_fasta_path)

        if aln_len == 0:
            print(f"[Skip] {filename} seems empty.")
            continue

        partition_filename = f"{gene_id}.nex"
        partition_path = os.path.join(gene_folder, partition_filename)
        create_partition_file(partition_path, aln_len)

        # D. Run IQ-TREE
        print(f"Processing: {gene_id}")
        print(f"  -> Length: {aln_len} bp (spaces removed)")
        print(f"  -> Running IQ-TREE...", end=" ", flush=True)

        cmd = [
            'iqtree',
            '-s', filename,             # Input (relative to cwd)
            '-sp', partition_filename,  # Partition (relative to cwd)
            '-m', 'MFP+MERGE',          # ModelFinder + Merge partitions
            '-nt', 'AUTO',              # Threads
            '-pre', gene_id,            # Prefix for output files
            '-bb', '1000'               # Ultrafast Bootstrap
        ]

        try:
            # We set cwd=gene_folder so all IQ-TREE output dumps into that folder
            subprocess.run(
                cmd,
                cwd=gene_folder,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                check=True
            )
            print("Done.")
        except subprocess.CalledProcessError as e:
            print("FAILED.")
            print(f"  [Error] IQ-TREE failed for {filename}")
            # print(e.stderr.decode())

    print("=" * 60)
    print(f"Pipeline complete. Data organized in: {output_root}/")

if __name__ == "__main__":

    # CONFIGURATION
    INPUT_ALIGNMENTS = "./combined_ortho_homologs/trimmed/"
    OUTPUT_DIR = "Tree_and_analyses"

    run_phylogeny_pipeline(INPUT_ALIGNMENTS, OUTPUT_DIR)
