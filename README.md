# 1. Download Orthologues from Ensembl
A script useful to download all orthologues as cds, of a Human Gene from the Compara database.
1_fetch_orthologues_2g.py

It uses ENSEMBL REST API.

My previous preferred way to download cds was (see here: https://pmc.ncbi.nlm.nih.gov/articles/PMC3855309/) using the EASER script, which used PyCogent, which is now not developed and has been modified to cogent3 and EnsemblLite (ensembl-tui or eti). This script can be modified for every use case in the rest.ensembl.org and download orthologues. 

```python3 fetch_orthologues.py```

The input file is given as a text file via sys_argv (see line 140: "input_arg = sys.argv[1]")

# 2. Align the downloaded sequences and Trim them

The script ```python3 2_6_align_and_trim.py```, is prepared to align the already downloaded sequences using MAFFT [^Katoh K, Misawa K, Kuma K, et al. MAFFT: a novel method for rapid multiple sequence alignment based on fast Fourier transform. Nucleic Acids Res 2002;30:3059–66.] and then trim the sequences considering them as codons. With the "relaxed" parameters of having half of the sequences with gaps or ambiguities. This same script is used in step 6 as well.

# 3. Finding Homologs of the downloaded "genes" from our transcriptomes

We can use the script ```python3 3_fetch_homologs_hmmer.py``` to 

* Make a HMMER (http://hmmer.org) profile using the alignments 
* Conduct HMMER-SEARCH against our CDS sequences (Transcriptome assembly) 
* output the homologue sequences into its own files

# 4. Create a "Homolog" fasta file join it with "ortholog" fasta file

If we are using multiple transcriptome assemblies for our study we can concatenate the "files" (not sequences) into one to create a single fasta file per gene. Which can be joined to the downloaded Orthologue fasta, to create a singe fasta that contains:

* Ensembl orthologs (from compara)
* Homologues (orthologues) from our assembly
We can use ```python3 4_conc_homologs.py``` and ```python3 5_merge_ensembl_and_homologs.py```

# 6. Align and trim the joined (orthologue from ensembl and our assembly) fasta file
# 7. Prepare ML phylogenies for the alignment
```python3 7_run_iqtree_pipeline.py```
we use:

* Ultrafast bootstrap 1000
* -m MFP+MERGE
* Partitioning by codon positions (1,2,3)
