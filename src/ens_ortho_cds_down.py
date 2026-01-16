#!/usr/bin/python3
# Importing the required libraries
import requests, sys
import time
import xml.etree.ElementTree as ET
import pandas as pd
import os, glob, shutil

#a script to download cds of a given human geneid using the ENSEMBL REST APIs
#my previous preferred way to download cds was (see here: https://pmc.ncbi.nlm.nih.gov/articles/PMC3855309/) using the EASER script, which used PyCogent,
#which is now not developed and has been modified to cogent3 and EnsemblLite (ensembl-tui or eti)
#this script can be modified for every use case in the rest.ensembl.org and download anything the user likes 
#usage python3 ensid.py name_of_the_text_file_with_ensemble_ids_one_id_per_line.txt
#definitions to download from restapi given below=============
def download_seqids_ensembl_restapi(ext_master):    
    server = "https://rest.ensembl.org"
    #ext_master
    ext=ext_master
    try:
        r = requests.get(server+ext, headers={ "Content-Type" : "text/xml"},timeout=60)
        file.write(r.text)
    except requests.exceptions.RequestException as e:
        print(e)
#downloading all ortholog information to xml file using the above function ===========================================================================
#getting transcript information as jsons
def download_transcript_ensembl_restapi(gid):
    server = "https://rest.ensembl.org"
    ext = '/overlap/translation/'+otd+'?'
    #print(server+ext)
    try:
        rs = requests.get(server+ext, headers={ "Content-Type" : "application/json"},timeout=60)
        file.write(rs.text)
    except requests.exceptions.RequestException as e:
        print(e)
#download cds sequences in fasta format
def download_seqs_ensembl_restapi(tsc):
    server = "https://rest.ensembl.org"
    ext = '/sequence/id/'+transcript_id+'?type=cds'
    #print(server+ext)
    try:
        rs = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"},timeout=60)
        file.write(rs.text)
    except requests.exceptions.RequestException as e:
        print(e)
#starting the script=========
masterdir=os.getcwd()#declaring the working directory
ensid_file = sys.argv[1]#the argument takes the file to read
with open(ensid_file, "r") as enfil:
    #ensid_list=[]#we open a list to store all the ensembl ids that we work with, intended to exit the for loop gracefully
    ensid_list = os.listdir()
    #ensid_list.append(dir_list)
    #print(ensid_list)
    for line in enfil:
        linestripped=line.splitlines()
        ensid=linestripped[0]
        #print(ensid)
        if ensid in ensid_list:
            pass
        elif ensid not in ensid_list:
            ensid_list.append(ensid)
           #print(ensid_list)
            os.makedirs(ensid)
            chdir=masterdir+'/'+ensid
            os.chdir(chdir)
            fn=ensid+'.xml'
        #print(fn)
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
                rows.append({"id": gid,"species": spp,"type": op,"pid": pid,"taxonomy_level": taxonomy_level})

            df = pd.DataFrame(rows, columns=cols)
# Writing dataframe to csv
            df.to_csv('output.csv')
            ndf = df[df['taxonomy_level']=='Euteleostomi']#selecting only fishes
#ndf = df[(df['taxonomy_level']=='Euteleostomi') & (df['taxonomy_level']=='Sarcopterygii')] # if we need latimeria can add & to get hag fish and petromyzon using 'vertebrata'
            ndf.to_csv('euteleostomi_orthologues.csv')
#Downloading json files for eveny pids to extract Parent===============================================================================
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
            path_json = paths+"/jsons"
            path_fasta = paths+"/fasta_cds"
            for file in files:#let's move all json file to the cerated 'jsons' directory
                if os.path.isfile(file):
                    shutil.move(file, path_json)
            for fasta_file in fastas:#Let's move all the fasta files to the created fasta_cds directory
                if os.path.isfile(fasta_file):
                    shutil.move(fasta_file, path_fasta)
            path_fasta_files = path_fasta+"/*.fasta"
            fasta_files = glob.glob(path_fasta_files)
            with open("concatenated_fasta.fasta", "wb") as outfile:#now we are going to concatenate the fasta files in the fasta_cds directory to a concatenated_fasta.fasta file
                for filename in fasta_files:
                    with open(filename,"rb") as infile:
                        shutil.copyfileobj(infile,outfile)
            os.chdir(masterdir)
#print("Orthologs for these human ensembl gene ids were downloaded: "+ensid_list)
