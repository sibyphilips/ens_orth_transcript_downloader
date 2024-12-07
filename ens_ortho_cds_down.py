#!/usr/bin/python3
# Importing the required libraries
import requests, sys
import time
import xml.etree.ElementTree as ET
import pandas as pd
import os, glob, shutil

#a script to download cds of a given human geneid using the ENSEMBL REST APIs

ensid=input("give a human geneid from ensembl to search: ")
def download_seqids_ensembl_restapi(ext_master):    
    server = "https://rest.ensembl.org"
    #ext_master
    ext=ext_master
    try:
        r = requests.get(server+ext, headers={ "Content-Type" : "text/xml"},timeout=60)
        file.write(r.text)
    except requests.exceptions.RequestException as e:
        print(e)

#downloading all orthologs using the above function ===========================================================================
fn=ensid+'.xml'
search_id = str(ensid)
ext_master = '/homology/id/human/'+search_id+'?'+'type=orthologues;format=condensed'#"/homology/id/human/ENSG00000157764?type=orthologues;format=condensed"
with open(fn, 'w') as file:
    download_seqids_ensembl_restapi(ext_master)
time.sleep(5) 
# Parsing the XML file=========================================================================================

cols = ["id", "species", "type", "pid", "taxonomy_level"]
rows = []
tree = ET.parse(fn)
root = tree.getroot()
speciesl = []
gidl = []
pidl = []
opl = []
for child in tree.iter('homologies'):
    spp = child.attrib['species']
    gid = child.attrib['id']
    pid = child.attrib['protein_id']
    op = child.attrib['type']
    taxonomy_level = child.attrib['taxonomy_level']
    rows.append({"id": gid,
                 "species": spp,
                 "type": op,
                 "pid": pid,
                 "taxonomy_level": taxonomy_level})

df = pd.DataFrame(rows, columns=cols)
# Writing dataframe to csv
df.to_csv('output.csv')
ndf = df[df['taxonomy_level']=='Euteleostomi']#selecting only fishes
#ndf = df[(df['taxonomy_level']=='Euteleostomi') & (df['taxonomy_level']=='Sarcoptrygii')] # if we need latimeria can add & to get hag fish and petromyzon using 'vertebrata'
ndf.to_csv('euteleostomi_orthologues.csv')
#Downloading json files for eveny pids to extract Parent===============================================================================
def download_transcript_ensembl_restapi(gid):
    server = "https://rest.ensembl.org"
    ext = '/overlap/translation/'+otd+'?'
    #print(server+ext)
    try:
        rs = requests.get(server+ext, headers={ "Content-Type" : "application/json"},timeout=60)
        file.write(rs.text)
    except requests.exceptions.RequestException as e:
        print(e)


def download_seqs_ensembl_restapi(tsc):
    server = "https://rest.ensembl.org"
    ext = '/sequence/id/'+transcript_id+'?type=cds'
    #print(server+ext)
    try:
        rs = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"},timeout=60)
        file.write(rs.text)
    except requests.exceptions.RequestException as e:
        print(e)

dfe = pd.read_csv('euteleostomi_orthologues.csv')
for gid in dfe['pid']:#only pid has parent transcript see https://www.biostars.org/p/105773/ which we collect below as tsc
    #print(gid)&spp in dfe['species']
    otd = str(gid)

    fln = otd+'.json'
    fln_fasta = otd+'.fasta'
    #server = "https://rest.ensembl.org"
    
    with open(fln, 'w') as file:
        download_transcript_ensembl_restapi(gid)
    #time.sleep(5)
    data_df = pd.read_json(fln) #reading the json file here
    tsc = data_df['Parent'].iloc[0]#extracting the transcript id here
    #print(tsc)
    transcript_id = str(tsc)
    fastaname = transcript_id+'.fasta'
    with open(fastaname,'w') as file:
        download_seqs_ensembl_restapi(transcript_id)


if not os.path.exists("jsons"):#creating a folder to store all json files
    os.makedirs("jsons")
if not os.path.exists("fasta_cds"):#creating a folder to store all fasta files
    os.makedirs("fasta_cds")
curr_dir = os.getcwd()#what is the current directory
files = glob.iglob(os.path.join(curr_dir, "*.json"))#make a list of all json files in this directory
fastas = glob.iglob(os.path.join(curr_dir, "*.fasta"))#make a list of all fasta files in this directory
paths = os.path.abspath(os.getcwd())#getting the path of the current directory to create a path for all json and fasta files to be copied in the next commands
#paths
#'/home/siby/workspace/works/ensembl_dl'
path_json = paths+"/jsons"
path_fasta = paths+"/fasta_cds"
#path_json
#'/home/siby/workspace/works/ensembl_dl/jsons'
for file in files:#let's move all json file to the cerated 'jsons' directory
    if os.path.isfile(file):
        shutil.move(file, path_json)

for fasta_file in fastas:#Let's move all the fasta files to the created fasta_cds directory
    if os.path.isfile(fasta_file):
        shutil.move(fasta_file, path_fasta)

path_fasta_files = path_fasta+"/*.fasta"
#path_fasta_files
#'/home/siby/workspace/works/ensembl_dl/fasta_cds/*.fasta'
fasta_files = glob.glob(path_fasta_files)
with open("concatenated_fasta.fasta", "wb") as outfile:#now we are going to concatenate the fasta files in the fasta_cds directory to a concatenated_fasta.fasta file
    for filename in fasta_files:
        with open(filename,"rb") as infile:
            shutil.copyfileobj(infile,outfile)
