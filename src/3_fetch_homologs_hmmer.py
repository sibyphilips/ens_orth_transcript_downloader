import os
import subprocess
import sys
import shutil

class HmmerPipeline:
    def __init__(self, assembly_path, evalue=1e-5, bin_path="/home/bin"):
        self.assembly_path = os.path.abspath(assembly_path)
        self.evalue = evalue
        self.bin_path = bin_path

        # Load the assembly sequences into memory once (for fast retrieval of best hits)
        print(f"Loading assembly sequences from {os.path.basename(assembly_path)}...")
        self.assembly_db = self._load_fasta_db(self.assembly_path)
        print(f"Loaded {len(self.assembly_db)} sequences.")

    def _load_fasta_db(self, fasta_path):
        """Helper: Reads the assembly fasta into a dictionary {header: sequence}"""
        db = {}
        header = None
        seq_lines = []
        try:
            with open(fasta_path, 'r') as f:
                for line in f:
                    if line.startswith(">"):
                        if header:
                            # Store previous sequence (take first word of header as key)
                            key = header.split()[0]
                            db[key] = "".join(seq_lines)
                        header = line.strip()[1:] # Remove '>'
                        seq_lines = []
                    else:
                        seq_lines.append(line.strip())
                # Store last sequence
                if header:
                    key = header.split()[0]
                    db[key] = "".join(seq_lines)
        except Exception as e:
            print(f"Error loading assembly: {e}")
            sys.exit(1)
        return db

    def _get_executable(self, tool_name):
        """Constructs path to executable in /home/bin"""
        # If the tool is directly in path, use it, otherwise check /home/bin
        path = os.path.join(self.bin_path, tool_name)
        if os.path.exists(path):
            return path
        else:
            # Fallback to system path if not in /home/bin
            return tool_name

    def build_hmm_profile(self, aligned_fasta, hmm_output):
        """Build HMM profile from aligned sequences"""
        hmmbuild_exe = self._get_executable('hmmbuild')
        try:
            # Assuming DNA since the search is nhmmer
            cmd = [hmmbuild_exe, '--dna', str(hmm_output), str(aligned_fasta)]

            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            return True, "HMM built"
        except subprocess.CalledProcessError as e:
            return False, f"hmmbuild error: {e.stderr}"
        except FileNotFoundError:
            return False, f"hmmbuild not found at {hmmbuild_exe}"

    def run_nhmmer(self, hmm_file, output_tbl):
        """Run nhmmer search (DNA sequences)"""
        nhmmer_exe = self._get_executable('nhmmer')
        try:
            cmd = [
                nhmmer_exe,
                '--tblout', str(output_tbl),
                '-E', str(self.evalue),
                str(hmm_file),
                str(self.assembly_path)
            ]

            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            return True, "Search completed"
        except subprocess.CalledProcessError as e:
            return False, f"nhmmer error: {e.stderr}"
        except FileNotFoundError:
            return False, f"nhmmer not found at {nhmmer_exe}"

    def extract_best_hit(self, tbl_file, output_fasta):
        """
        Parses the --tblout file to find the best hit (lowest E-value).
        Extracts that sequence from self.assembly_db.
        """
        best_hit_name = None
        best_evalue = float('inf')

        try:
            with open(tbl_file, 'r') as f:
                for line in f:
                    if line.startswith("#"): continue

                    parts = line.split()
                    if len(parts) < 13: continue

                    # nhmmer tblout format:
                    # target_name (0) ... evalue (12) ...
                    target_name = parts[0]
                    try:
                        evalue = float(parts[12])
                    except ValueError:
                        continue

                    # We want the lowest E-value (Best Hit)
                    if evalue < best_evalue:
                        best_evalue = evalue
                        best_hit_name = target_name

            if best_hit_name:
                if best_hit_name in self.assembly_db:
                    sequence = self.assembly_db[best_hit_name]
                    with open(output_fasta, 'w') as outf:
                        outf.write(f">{best_hit_name} [Best Hit E={best_evalue}]\n{sequence}\n")
                    return True, f"Found hit: {best_hit_name} (E={best_evalue})"
                else:
                    return False, f"Hit found in table ({best_hit_name}) but sequence missing in FASTA DB."
            else:
                return False, "No significant hits found."

        except Exception as e:
            return False, f"Error parsing table: {e}"

# ==========================================
# MAIN WORKFLOW
# ==========================================
def main():
    # --- CONFIGURATION ---
    # Path to the folder containing your aligned/trimmed fasta files
    INPUT_TRIMMED_DIR = "./Downloads/trimmed/"

    # Path to your Transcriptome Assembly Fasta
    # CHANGE THIS to the actual path of your assembly file
    ASSEMBLY_FILE = "/run/media/siby/TOSHIBA EXT/Transcriptome_Bini/3.Assembly/SD_trinity.Trinity.cdhit.fasta"

    # =====================

    if not os.path.exists(ASSEMBLY_FILE):
        print(f"Error: Assembly file not found at {ASSEMBLY_FILE}")
        return

    # 1. Setup Directories
    hmm_profile_dir = "hmm_profiles"
    if not os.path.exists(hmm_profile_dir):
        os.makedirs(hmm_profile_dir)

    # Create output folder named after the assembly file (without extension)
    assembly_name = os.path.splitext(os.path.basename(ASSEMBLY_FILE))[0]
    results_dir = f"{assembly_name}_hits"
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    print(f"Output directories:\n  Profiles: {hmm_profile_dir}\n  Results:  {results_dir}\n")

    # 2. Initialize Pipeline
    pipeline = HmmerPipeline(ASSEMBLY_FILE, evalue=1e-5)

    # 3. Process Files
    fasta_files = [f for f in os.listdir(INPUT_TRIMMED_DIR) if f.endswith(".fasta") or f.endswith(".fa")]

    if not fasta_files:
        print(f"No fasta files found in {INPUT_TRIMMED_DIR}")
        return

    print(f"Processing {len(fasta_files)} alignments...")
    print("-" * 60)

    for filename in fasta_files:
        base_name = os.path.splitext(filename)[0]
        input_path = os.path.join(INPUT_TRIMMED_DIR, filename)

        print(f"Gene: {base_name}")

        # --- Step A: Build HMM Profile ---
        hmm_path = os.path.join(hmm_profile_dir, f"{base_name}.hmm")
        ok_build, msg_build = pipeline.build_hmm_profile(input_path, hmm_path)

        if not ok_build:
            print(f"  [X] HMM Build Failed: {msg_build}")
            continue

        # --- Step B: Run Search (nhmmer) ---
        tbl_path = os.path.join(results_dir, f"{base_name}_hits.tbl")
        ok_search, msg_search = pipeline.run_nhmmer(hmm_path, tbl_path)

        if not ok_search:
            print(f"  [X] Search Failed: {msg_search}")
            continue

        # --- Step C: Extract Best Hit ---
        hit_fasta_path = os.path.join(results_dir, f"{base_name}_best_hit.fasta")
        ok_extract, msg_extract = pipeline.extract_best_hit(tbl_path, hit_fasta_path)

        if ok_extract:
            print(f"  [V] {msg_extract}")
        else:
            print(f"  [-] {msg_extract}")

    print("-" * 60)
    print("Pipeline complete.")

if __name__ == "__main__":
    main()
