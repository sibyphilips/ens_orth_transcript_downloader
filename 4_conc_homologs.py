import os
import sys

def generate_output_filename(original_filename):
    """
    Creates a new filename based on the first two words of the original file.
    Example:
      Input:  "ABHD11_ENSG00000106077_fishes.fasta"
      Output: "ABHD11_ENSG00000106077.fasta"
    """
    base_name = os.path.splitext(original_filename)[0]
    normalized_name = base_name.replace(" ", "_")
    parts = normalized_name.split("_")

    if len(parts) >= 2:
        new_name = f"{parts[0]}_{parts[1]}.fasta"
    else:
        new_name = f"{parts[0]}.fasta"

    return new_name

def get_folder_suffix(folder_path):
    """
    Extracts the first word of the folder name (separated by underscores).
    Example: "./human_homologs" -> suffix "human"
    """
    # Normalize path to handle trailing slashes
    clean_path = os.path.normpath(folder_path)
    # Get the directory name (e.g., "human_homologs")
    folder_name = os.path.basename(clean_path)
    # Get the first word before any underscore
    first_word = folder_name.split('_')[0]
    return first_word

def concatenate_homologs(folder_list, output_dir):
    """
    Concatenates FASTA files, cleans headers, and adds folder-specific suffixes.
    """
    # 1. Validate Input
    if not folder_list:
        print("Error: No input folders provided.")
        return

    for f in folder_list:
        if not os.path.exists(f):
            print(f"Error: Input folder does not exist: {f}")
            return

    # 2. Create Output Directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 3. Get Master List of Files from First Folder
    reference_folder = folder_list[0]
    reference_files = [
        f for f in os.listdir(reference_folder)
        if f.endswith(('.fasta', '.fa', '.fas'))
    ]

    print(f"Found {len(reference_files)} files in reference folder: {os.path.basename(reference_folder)}")
    print("-" * 60)

    # 4. Process Each File
    count_processed = 0

    for filename in reference_files:
        output_filename = generate_output_filename(filename)
        output_path = os.path.join(output_dir, output_filename)

        found_count = 0

        with open(output_path, 'w') as outfile:

            # Iterate through EVERY folder in the list
            for folder in folder_list:
                source_file_path = os.path.join(folder, filename)

                if os.path.exists(source_file_path):
                    found_count += 1

                    # Determine suffix for this specific folder
                    # e.g., "human" from "human_homologs"
                    suffix = get_folder_suffix(folder)

                    try:
                        with open(source_file_path, 'r') as infile:
                            for line in infile:
                                if line.startswith(">"):
                                    # HEADER MODIFICATION LOGIC

                                    # 1. Strip whitespace
                                    clean_line = line.strip()

                                    # 2. Remove square brackets info (Take everything before first space or bracket)
                                    # Example: ">TRINITY... [Best Hit]" -> ">TRINITY..."
                                    # Using split()[0] effectively grabs just the ID part
                                    header_id = clean_line.split()[0]

                                    # 3. Construct new header with suffix
                                    # Output: >TRINITY..._human
                                    new_header = f"{header_id}_{suffix}\n"
                                    outfile.write(new_header)

                                else:
                                    # SEQUENCE LINE
                                    outfile.write(line)

                            # Ensure newline at end of file to prevent merging
                            outfile.write('\n')

                    except Exception as e:
                        print(f"  [Error] Reading {source_file_path}: {e}")
                else:
                    print(f"  [Warning] File {filename} missing in {os.path.basename(folder)}")

        count_processed += 1

    print("-" * 60)
    print(f"Successfully processed {count_processed} gene sets.")
    print(f"Combined files saved to: {output_dir}")

# ==========================================
# CONFIGURATION
# ==========================================
if __name__ == "__main__":

    # --- 1. LIST YOUR INPUT FOLDERS HERE ---
    # Ensure these names are descriptive (e.g., "human_data", "fish_homologs")
    # because the FIRST WORD will be used as the suffix.
    INPUT_FOLDERS = [
        "./DF_trinity.Trinity.cdhit_hits",      # Suffix will be "_human"
        "./HF_trinity.Trinity.cdhit_hits",      # Suffix will be "_mouse"
        "./SD_trinity.Trinity.cdhit_hits",  # Suffix will be "_zebrafish"
    ]

    # --- 2. OUTPUT FOLDER NAME ---
    OUTPUT_FOLDER = "combined_homologs"

    concatenate_homologs(INPUT_FOLDERS, OUTPUT_FOLDER)
