# ens_orth_transcript_downloader
A script useful to download all orthologues as cds, of a Human Gene from the Compara database.

# Uses ENSEMBL REST API
# Downloads cds of orthologues in compara database
This circumvents the step of finding the orthologues of all the genes that you are interested using some external programme. The ensembl database has the orthologues already calculated and stored in the compara database.
This should be useful for people looking to work on comparative genomics projects.

# steps
- A user is asked to provide a Human Geneid for the corresponding orthologues that the user wants.
- All orthologues (one to one, one to many, many to many etc.) are downloaded from the ENSEMBL COMPARA database
- The resulting *.csv file is parsed to get the protein IDs and then search for the TRANSCRIPTS
- The transcript IDs are used to download the corresponding CDS as separate fasta files
- JSON files are downloaded and shifted to a json file folder
- Fasta files are downloaded per CDS and shifted to a fasta file folder
- A Single Fasta File containing all the downloaded CDS is produced in the current working directory
- An XML file is first prepared containing all the orthologues of the provided geneid
- This xml file is parsed as csv files (here since I am interested on fishes I have made a separate csv file containing only the 'euteleostomi' orthologues)

# Usage
```$python3 /path/to/the/script/this_script.py```

```$give a human geneid from ensembl to search:#provide the ENSG--- here ```
