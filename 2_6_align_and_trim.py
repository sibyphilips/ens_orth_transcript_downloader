import subprocess
import os
import shutil
import re

# ==========================================
# 1. HELPER: SEQUENCE MAPPER CLASS
# ==========================================
class SequenceMapper:
    """
    Handles renaming FASTA headers to safe, short IDs (e.g., >Seq_001)
    to appease GBlocks, and restoring original headers afterwards.
    """
    def __init__(self):
        self.mapping = {}  # Maps "Seq_001" -> "Original_Long_Header..."

    def create_temp_safe_fasta(self, input_path, temp_path):
        """
        Reads input_path, writes temp_path with short headers (>Seq_X).
        Stores the mapping.
        """
        self.mapping = {}
        count = 0
        try:
            with open(input_path, 'r') as f_in, open(temp_path, 'w') as f_out:
                for line in f_in:
                    if line.strip().startswith(">"):
                        original_header = line.strip()[1:] # Remove '>'

                        # Create safe ID
                        count += 1
                        safe_id = f"Seq_{count}"

                        # Store map
                        self.mapping[safe_id] = original_header

                        # Write safe header
                        f_out.write(f">{safe_id}\n")
                    else:
                        f_out.write(line)
            return True, count
        except Exception as e:
            return False, f"Error creating temp file: {e}"

    def restore_original_headers(self, safe_file_path, final_output_path):
        """
        Reads the GBlocks output (safe_file_path) which has >Seq_X headers,
        and writes final_output_path with Original headers.
        """
        try:
            with open(safe_file_path, 'r') as f_in, open(final_output_path, 'w') as f_out:
                for line in f_in:
                    if line.strip().startswith(">"):
                        # Extract the ID (e.g., Seq_1)
                        # GBlocks sometimes adds spaces or info after the ID, so we split
                        safe_id_candidate = line.strip()[1:].split()[0]

                        if safe_id_candidate in self.mapping:
                            original = self.mapping[safe_id_candidate]
                            f_out.write(f">{original}\n")
                        else:
                            f_out.write(line)
                    else:
                        f_out.write(line)
            return True, "Success"
        except Exception as e:
            return False, f"Error restoring headers: {e}"

# ==========================================
# 2. CORE FUNCTIONS
# ==========================================

def run_mafft(input_file, output_file):
    """ Runs MAFFT alignment using absolute paths. """
    input_abs = os.path.abspath(input_file)
    output_abs = os.path.abspath(output_file)

    cmd = ['mafft', '--auto', input_abs]

    try:
        with open(output_abs, 'w') as outf:
            subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, text=True, check=True)
        return True, "Success"
    except subprocess.CalledProcessError as e:
        return False, f"MAFFT Error: {e.stderr}"
    except FileNotFoundError:
        return False, "MAFFT not found."

def run_gblocks_safely(aligned_file, output_folder, mol_type='c'):
    """
    1. Creates a temp file with short names.
    2. Runs GBlocks.
    3. Restores original names into *_aln_tr.fasta.
    4. Moves HTML report to trimmed folder as *_aln_tr.html.
    """
    aligned_abs = os.path.abspath(aligned_file)

    # Define temp file paths (in same directory as aligned file)
    temp_safe_input = aligned_abs + ".temp_safe"

    # Initialize Mapper
    mapper = SequenceMapper()

    # A. Create Safe Temp File
    ok, result = mapper.create_temp_safe_fasta(aligned_abs, temp_safe_input)
    if not ok:
        return False, result

    num_seqs = result
    if num_seqs == 0:
        return False, "Aligned file contains 0 sequences."

    print(f"    -> GBlocks Input: {os.path.basename(aligned_file)} ({num_seqs} sequences)")

    # B. Calculate Parameters
    b1_val = int(round((51 / 100) * num_seqs))
    b2_val = int(round((85 / 100) * num_seqs))
    if b1_val > b2_val: b2_val = b1_val

    # C. Run GBlocks
    # We use extension ".gb"
    result_suffix = ".gb"
    cmd = [
        'Gblocks',
        temp_safe_input,
        f'-t={mol_type}',
        f'-b1={b1_val}',
        f'-b2={b2_val}',
        '-b3=5',
        '-b4=10',
        '-b5=h',
        f'-e={result_suffix}'
    ]

    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # --- FIX: CORRECT FILENAME PREDICTION ---
        # GBlocks output will be: inputfilename + suffix
        # GBlocks HTML will be:   inputfilename + suffix + .htm

        temp_gblocks_output = temp_safe_input + result_suffix          # e.g., file.temp_safe.gb
        temp_gblocks_html   = temp_safe_input + result_suffix + ".htm" # e.g., file.temp_safe.gb.htm

        if os.path.exists(temp_gblocks_output):
            # --- NAMING LOGIC ---
            base_name = os.path.splitext(os.path.basename(aligned_file))[0]
            if base_name.endswith("_aligned"):
                clean_name = base_name[:-8]
            else:
                clean_name = base_name

            final_fasta_name = f"{clean_name}_aln_tr.fasta"
            final_html_name = f"{clean_name}_aln_tr.html"

            final_fasta_path = os.path.join(output_folder, final_fasta_name)
            final_html_path = os.path.join(output_folder, final_html_name)

            # --- D. Restore Headers (Save FASTA) ---
            restore_ok, restore_msg = mapper.restore_original_headers(temp_gblocks_output, final_fasta_path)

            # --- E. Move HTML File ---
            if os.path.exists(temp_gblocks_html):
                shutil.move(temp_gblocks_html, final_html_path)
            else:
                # Debugging fallback: check if it exists without the .gb middle part?
                # (Some versions vary, but .gb.htm is standard behavior for -e=.gb)
                alt_html = temp_safe_input + ".htm"
                if os.path.exists(alt_html):
                    shutil.move(alt_html, final_html_path)

            # --- F. Cleanup Temp Files ---
            # Remove the safe input
            if os.path.exists(temp_safe_input):
                os.remove(temp_safe_input)
            # Remove the raw GBlocks output (we have the restored version now)
            if os.path.exists(temp_gblocks_output):
                os.remove(temp_gblocks_output)
            # Remove any leftover HTML if move failed
            if os.path.exists(temp_gblocks_html):
                os.remove(temp_gblocks_html)

            if restore_ok:
                return True, "Success"
            else:
                return False, f"Error restoring headers: {restore_msg}"
        else:
            # Cleanup input temp if failed
            if os.path.exists(temp_safe_input): os.remove(temp_safe_input)

            return False, f"GBlocks produced no output.\nSTDOUT: {proc.stdout}"

    except Exception as e:
        return False, f"Execution Error: {e}"

# ==========================================
# 3. MAIN PIPELINE
# ==========================================

def batch_process():
    # SETTINGS
    INPUT_FOLDER = "."
    MOLECULE_TYPE = "c"  # c = Codons

    # Setup Folders
    aligned_dir = os.path.join(INPUT_FOLDER, "aligned")
    trimmed_dir = os.path.join(INPUT_FOLDER, "trimmed")

    for d in [aligned_dir, trimmed_dir]:
        if not os.path.exists(d): os.makedirs(d)

    # Find Files
    all_files = os.listdir(INPUT_FOLDER)
    fasta_files = []
    for f in all_files:
        if f.lower().endswith(".fasta") and "_aligned" not in f and "_aln_tr" not in f:
            fasta_files.append(f)

    if not fasta_files:
        print("No .fasta files found.")
        return

    print(f"Found {len(fasta_files)} FASTA files. Processing...")
    print("-" * 60)

    for filename in fasta_files:
        input_path = os.path.join(INPUT_FOLDER, filename)
        base_name = os.path.splitext(filename)[0]

        # 1. ALIGN
        aligned_filename = f"{base_name}_aligned.fasta"
        aligned_path = os.path.join(aligned_dir, aligned_filename)

        print(f"Processing: {filename}")
        print(f"  1. Aligning...", end=" ", flush=True)

        mafft_ok, mafft_msg = run_mafft(input_path, aligned_path)

        if not mafft_ok:
            print(f"FAILED. {mafft_msg}")
            continue
        print("Done.")

        # 2. TRIM
        print(f"  2. Trimming...", end="\n", flush=True)

        gb_ok, gb_msg = run_gblocks_safely(aligned_path, trimmed_dir, mol_type=MOLECULE_TYPE)

        if gb_ok:
            print("     -> Done. Saved to 'trimmed/' (*_aln_tr.fasta & *.html)")
        else:
            print(f"     -> FAILED. {gb_msg}")

        print("-" * 60)

if __name__ == "__main__":
    batch_process()
