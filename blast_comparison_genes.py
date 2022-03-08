#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 11:58:22 2021

@author: trucdo

This script takes a list of genes of interest from the B. cepacia ATCC 25416 genome 
and extracts the BLASTp results for those genes against the other Bcc strains. 
The BLASTp results are shown as a heat map.

input: .txt files of BLASTp results against each strain (in -outfmt "7" command line option)
input: .txt file containing newline-delimited list of subject organisms used in BLAST
input: .txt file of protein accession nos. to compare across BLASTp result files
output: .fasta file of protein homologs for a query protein of interest across all subject organisms
define: query_organism; subject_organism_filename; query_proteins_filename; io_file_path
define: io_subfolder; blast_file_name_format; qcovhsp_cutoff; score_type
note: cross-check column header of BLASTp results file with blastp_record class definition

"""

#%% Define path, file names, and other parameters

# name of query organism that was BLASTed
query_organism = "Burkholderia_cepacia_ATCC_25416"

# name of .txt file containing list of subject organisms (databases) against which the query organism was individually BLASTed
subject_organism_filename = "subject_organism_names.txt"

# name of .txt file containing query organism genes to search
query_proteins_filename = "220107_Predicted_TBDTs.txt"

# main directory of input and output files
io_file_path = "path_to_main_directory"

# subpath to input files and where to print output files
io_subfolder =  io_file_path + "/" + query_organism + "/comparative_blast_outfmt7"

# BLASTp results files are named using the following format
# Code will replace <query_organism> and <subject_organism> with query organism and subject organism names, respectively
blast_file_name_format = "<query_organism>-<subject_organism>_blastdb.txt"

# indicate cutoff for percent coverage (qcovhsp) of alignment HSP
# cutoff will remove under-coveraged hits prior to heat map visualization to avoid skewing results
qcovhsp_cutoff = 50

# indicate parameter (str) to use for calculating protein similarity
# "pident" for percent identity, "ppos" for chemically-similar residues, "pidentpos" for average of pident and ppos
score_type = "ppos"

#%% Import modules 

# allows user to interact with OS-dependent functionality
import os
# for regular expression searches
import re
# import SeqIO subpackage from Bio package to handle FASTA file
from Bio import SeqIO
# to work with pandas dataframe
import pandas as pd
# import seaborn and matplotlib packages for plotting data
import seaborn
import matplotlib.pyplot as plt
# remember to import appropriate machine-learning packages in their own sections

#%% Class definitions

class fasta_record(object):
# define class to process FASTA-formatted protein record: each instance corresponds to one protein
# all attributes except length are strings
# input: SeqRecord object for one protein
    
    # constructor retrieves properties of protein record
    def __init__(self, protein):
        
        # protein accession number
        self.qacc = protein.id
        # full description of protein including accession number and annotated Genbank function
        self.description = protein.description
        # protein sequence
        self.seq = str(protein.seq)
        # length of protein
        self.length = len(protein.seq)



class blastp_record(object):
# define blastp_record class: each instance is the blastp result(s) for a query protein
# note that a query protein may have no hits or multiple alignments and multiple HSPs per alignment
# remember to import re module for blastp_record class definition to work
# argument: Blastp record for one protein in -outfmt 7 (str)
# note that -outfmt 7 produces an empty record labeled with "0 hits found" for queries with no hits
# all attributes are data-type str

    # constructor retrieves Blastp results and defines attributes
    def __init__(self, blast_record):
        
        # split the record for current protein using # delimiter separating the sections
        # if query protein had hits, then list is [Query, Database, Fields, num hits found]
        # if query protein had no hits, then list is [Query, Database, num hits found]
        # remember to remove newline character at the end of each item in the list
        blast_record_items = blast_record.split("# ")
        # remove the first item in the list, which is a blank item after splitting with delimiter
        blast_record_items.pop(0)
        
        # qacc and description of query protein (index 0 in list)
        query = blast_record_items[0].rstrip("\n")
        # split using ": " delimiter to retrieve only value (index 1) and not descriptor (index 0)
        self.query = query.split(": ")[1]
        
        # path of Blast database used
        dbpath = blast_record_items[1].rstrip("\n")
        self.db_path = dbpath.split(": ")[1]
        
        # name of Blast database used
        db_match_obj = re.search(r"/[A-Za-z0-9_]+_blastdb", self.db_path)
        blastdb = db_match_obj.group()
        self.blastdb = blastdb.rstrip("_blastdb").strip("/")
        
        # item at index 2 in list is either...
        # "Fields: header,..." indicating types of data in each tab-delimited column for queries with hits
        # "0 hits found" for queries with no hits
        item_index_2 = blast_record_items[2].rstrip("\n")
        
        # to pass to hsp_dict attribute: stores data for top HSP of top alignment
        # keys are renamed headers and values are data in corresponding column
        hsp_dict = {}
        
        # for query with no Blast hits
        if item_index_2 == "0 hits found":
            
            self.fields = "na"
            self.renamed_fields = "na"
            self.num_hits = "0"
            
            # split qacc and description of query protein using space-character delimiter
            # qacc is index 0 and description is index 1 in list
            qacc = self.query.split(" ")[0]
            
            # assign default values for all data columns because query had no hit
            hsp_dict["qacc"] = qacc
            hsp_dict["saccver"] = "no hit"
            hsp_dict["stitle"] = "no hit"
            hsp_dict["evalue"] = "10000.0" # arbitrary high number, convert to float as needed
            hsp_dict["bitscore"] = "0.0"
            hsp_dict["pident"] = "0.0"
            hsp_dict["ppos"] = "0.0"
            hsp_dict["qcov"] = "0.0"
            hsp_dict["qcovhsp"] = "0.0"
            hsp_dict["align_len"] = "0" # convert to int as needed
            hsp_dict["gap_opens"] = "na"
            hsp_dict["num_gaps"] = "na"
            hsp_dict["q_start"] = "na"
            hsp_dict["q_end"] = "na"
            hsp_dict["s_start"] = "na"
            hsp_dict["s_end"] = "na"
            hsp_dict["pidentpos"] = "0.0" # average of pident and ppos scores
        
        # for query with one or more Blast hits                  
        else:
            
            # split using ": " delimiter to retrieve only field names (index 1) and not descriptor (index 0)          
            self.fields = item_index_2.split(": ")[1]
            
            # rename each of the fields to remove periods and spaces and for ease of calling
            # if more data columns are added to Blastp files in the future, this makes it easy to add these data
            renamed_fields = self.fields
    
            if re.search(r"query acc\.", renamed_fields):
                renamed_fields = renamed_fields.replace("query acc.", "qacc")
            if re.search(r"subject acc\.ver", renamed_fields):
                renamed_fields = renamed_fields.replace("subject acc.ver", "saccver")
            if re.search(r"subject title", renamed_fields):
                renamed_fields = renamed_fields.replace("subject title", "stitle")
            if re.search(r"bit score", renamed_fields):
                renamed_fields = renamed_fields.replace("bit score", "bitscore")         
            if re.search(r"\% identity", renamed_fields):
                renamed_fields = renamed_fields.replace("% identity", "pident") 
            if re.search(r"\% positives", renamed_fields):
                renamed_fields = renamed_fields.replace("% positives", "ppos")   
            if re.search(r"\% query coverage per subject", renamed_fields):
                renamed_fields = renamed_fields.replace("% query coverage per subject", "qcov")  
            if re.search(r"\% query coverage per hsp", renamed_fields):
                renamed_fields = renamed_fields.replace("% query coverage per hsp", "qcovhsp")   
            if re.search(r"alignment length", renamed_fields):
                renamed_fields = renamed_fields.replace("alignment length", "align_len")
            if re.search(r"gap opens", renamed_fields):
                renamed_fields = renamed_fields.replace("gap opens", "gap_opens")   
            if re.search(r"gaps", renamed_fields):
                renamed_fields = renamed_fields.replace("gaps", "num_gaps")
            if re.search(r"q\. start", renamed_fields):
                renamed_fields = renamed_fields.replace("q. start", "q_start")
            if re.search(r"q\. end", renamed_fields):
                renamed_fields = renamed_fields.replace("q. end", "q_end")     
            if re.search(r"s\. start", renamed_fields):
                renamed_fields = renamed_fields.replace("s. start", "s_start")
            if re.search(r"s\. end", renamed_fields):
                renamed_fields = renamed_fields.replace("s. end", "s_end")   
            
            # split field names using comma delimiter
            fields_list = renamed_fields.split(", ")
            
            # tab-delimited renamed field names
            self.renamed_fields = "\t".join(fields_list)
            
            # item at index 3 of list for proteins with one or more hits is the section containing...
            # ...number of hits found and each subsequent line is the result of a hit
            # remove newline after last hit and then split hits using interior newline characters
            hits_list = blast_record_items[3].rstrip("\n").split("\n")
            # remove the first item in hits_list, which is number of hits found
            hits_list.pop(0)
            # now length of hits_list is the number of hits per protein query (combines HSPs for all alignments)
            self.num_hits = len(hits_list)
            
            # position of evalue column
            evalue_index = fields_list.index("evalue")
            # set an arbitrary high evalue cutoff to compare all hits for most significant hit per query
            e_bound = 10000.0 
            # iterate through each hit
            for hit in hits_list:
                # split the tab-delimited data for current hit in loop
                hit_split = hit.split("\t")
                # re-assign top-scoring HSP if its evalue is lower than the current evalue cutoff
                if float(hit_split[evalue_index]) < e_bound:
                    top_hsp_split = hit_split
                    e_bound = float(hit_split[evalue_index])
                
            # store field name (key) and corresponding data (value) from top HSP of top alignment in a dict
            counter = 0
            for field in fields_list:
                hsp_dict[field] = top_hsp_split[counter]
                counter = counter + 1
                
            # calculate and store the average of pident and ppos scores
            pidentpos = (float(hsp_dict["pident"]) + float(hsp_dict["ppos"])) / 2
            hsp_dict["pidentpos"] = str(pidentpos)
        
        # attribute which references a dict that stores data from the top HSP of top alignment
        self.hsp_dict = hsp_dict

#%% Function definitions

def process_blast_file(file_name, file_path):
# define a function to convert an -outfmt 7 blastp results file into a dict
# remember to import os module and define blastp_record class for process_blast_file function to work
# arguments: name of blast results file (str), path to blastp results file (str)
# return: dict object; key is query accession no. and value is a blastp_record object
    
    # navigate to folder containing blastp results
    os.chdir(file_path)
    blast_results_file_object = open(file_name)
    blast_results_file = blast_results_file_object.read()
    blast_results_file_object.close()
    
    # separate results found for each query protein (a query may have multiple alignments and HSPs)
    # use the delimiter "# BLASTP 2.10.0" which separates the record for each query
    blast_results_split = blast_results_file.split("# BLASTP 2.10.0+\n")
    
    # remove the first item in the list, which is a blank item after splitting with delimiter
    blast_results_split.pop(0)
    
    # dict-type variable to store and return parsed blast results
    process_blast_file_dict = {}
    
    # iterate through list of blastp results for each protein query
    # convert blastp results of each query to a blastp_record object
    # store the blastp_record objects as values in a dict with corresponding qacc as the keys
    for record in blast_results_split:
        temp_blast_obj = blastp_record(record)
        process_blast_file_dict[temp_blast_obj.hsp_dict["qacc"]] = temp_blast_obj
        
    return process_blast_file_dict
  

    
def process_strains(subj_org_filename, query_org, file_path, file_name_format):
# define a function to convert blastp results file into a dict of blastp_record objects for each strain
# remember to define process_blast_file function for process_strains function to work
# arguments: name of file containing newline-delimited list of subject organisms (str)
# arguments: name of query organism (str), Blast results file path (str), Blast results file name format (str)
# return: tuple of strains_list, strain_dict
    
    # navigate to folder containing blastp results
    os.chdir(file_path)
    strains_file_object = open(subj_org_filename)
    strains_file = strains_file_object.read()
    strains_file_object.close()
    
    # store each strain name as item in a list
    strains_list = strains_file.split("\n")
    
    # dict stores blastp results against current subject strain
    # key is strain name and value is inner dict of blastp_record objects (key: qacc; value: blastp_record objects)
    strains_dict = {}
    
    # iterate through each strain in the list and turn blastp results into 
    # an outer dict of strains:process_blast_file_dict within which is an
    # inner dict of of qacc:blastp_obj
    # replace the file_name_format variable with the appropriate organism names to open the correct files
    for strain in strains_list:
        current_blast_file = file_name_format.replace("<query_organism>", query_org)
        current_blast_file = current_blast_file.replace("<subject_organism>", strain)
        strains_dict[strain] = process_blast_file(current_blast_file, file_path)
    
    # return tuple of strains list, strains dict
    return (strains_list, strains_dict)



def search_blast(strains_dict, query_protein_filename, file_path):
# function to search qacc from protein list against all blastp results file and retrieve matching record
# arguments: dict of strains and their blastp_record objects (dict)
# arguments: name of file containing newline-delimited list of proteins to search (str)
# arguments: results file path (str)
# return: POI_blast_dict

    # navigate to folder containing blastp results
    os.chdir(file_path)
    proteins_file_object = open(query_protein_filename)
    proteins_file = proteins_file_object.read()
    proteins_file_object.close()
    
    # store each protein accession no. as item in a list
    proteins_list = proteins_file.split("\n")
    
    # multi-layered dict object stores blastp results against all subject strain
    # first layer: POI_dict[qacc] = strains_dict
    # second layer: strains_dict[strain_name] = blastp_record object for qacc for this strain
    POI_blast_dict = {}
    
    # iterate through each query protein of interest and retrieve blastp result using its qacc
    # against all subject organisms
    for protein in proteins_list:
        
        # temporary dict to store blastp of current protein qacc across all subject strains
        strain_subdict = {}
        
        for strain, qacc_dict in strains_dict.items():
            
            # retrieve and store the blastp result object matching current qacc for current strain
            qacc_match_blastp_obj = qacc_dict.get(protein, "no hit")
            
            # if query protein has no match, print to screen protein qacc and skip current protein 
            if qacc_match_blastp_obj == "no hit":
                print(strain)
                print(protein)
                break
            
            strain_subdict[strain] = qacc_match_blastp_obj
        
        if qacc_match_blastp_obj != "no hit":
            # store blastp result of current protein of interest against all subject strain
            POI_blast_dict[protein] = strain_subdict
            
    return POI_blast_dict


         
def dict_to_df(POI_blast_dict, score_type="pident", coverage_cutoff=0.0):
# define a function to convert score for all protein queries from all strains into a 2-D array
# remember to import pandas module for dict_to_df function to work
# arguments: compiled blastp results (dict; output of search_blast function)
# arguments: blastp score parameter (default pident, ppos, or average) to use for comparing hits (str)
# arguments: percent coverage qcovhsp cutoff for proteins with sufficient alignment (int)
# output: a 2-D Pandas DataFrame; each row is a protein query and each column is a different strain
# note that identity score values are of data-type float

    # obtain sorted list of strain names
    all_qacc_list = sorted(POI_blast_dict.keys())
    strain_list = []
    # iterate through first value (first protein query) in the POI_blast_dict    
    for strain, blastp_record in POI_blast_dict[all_qacc_list[0]].items():
        strain_list.append(strain)
               
    # create an empty 2-D Pandas DataFrame to later store scores of blastp hits
    # rows = qacc; columns = strain
    blast_df = pd.DataFrame(index = all_qacc_list, columns = strain_list, dtype = "float")

    # iterate through the dataframe query accession no. and search the compile_blast_dict
    for qacc in blast_df.index:
        
        # obtain all blastp_record objects for current qacc in loop
        strains_dict = POI_blast_dict[qacc]
        
        # iterate through each strain in dataframe
        for strain_name in blast_df.columns:
           
            # obtain blastp_record object for current strain in loop
            blast_obj = strains_dict[strain_name]
            # store identity score in cell located at row "qacc" and column "strain_name"
            #blast_df.at[qacc, strain_name] = float(blast_obj.hsp_dict[score_type])
            blast_df.at[qacc, strain_name] = float(blast_obj.hsp_dict[score_type])*(float(blast_obj.hsp_dict["qcovhsp"]) / 100)
    
    # reorder blast_df so that the genes are listed from top to bottom in increasing gene numbers  
    blast_df = blast_df.sort_index(axis=0, ascending = True)
        
    # reorder blast_df so that susceptible and non-susceptible strains are grouped separately
    column_names = ["Burkholderia_vietnamiensis_AU3578", \
                    "Burkholderia_vietnamiensis_AU21726", \
                        "Burkholderia_ambifaria_AMMD", \
                            "Burkholderia_multivorans_AU15814", \
                                "Burkholderia_dolosa_PC543", \
                                    "Burkholderia_dolosa_AU0158", \
                                        "Burkholderia_multivorans_ATCC_17616", \
                                            "Burkholderia_cenocepacia_AU0756", \
                                                "Burkholderia_cenocepacia_AU24362", \
                                                    "Burkholderia_cenocepacia_J2315", \
                                                        "Burkholderia_multivorans_AU30438"]
        
    blast_df = blast_df.reindex(columns = column_names)
  
    return blast_df


        
def fasta_to_dict(strains_list, file_path):
# define a function to convert FASTA files into a dict of fasta_record objects for each strain
# make sure to import SeqIO module from Bio package and os module and define fasta_record class
# arguments: list of organism names (list of str), main file directory (str)
# return: strains_fasta_dict
   
    # create a multi-layered dict to store fasta_record objects for all strains
    # outer dict: strain-inner dict
    # inner dict: protein accession no-fasta_record object
    strains_fasta_dict = {}
    
    for strain in strains_list:
        
        # navigate to folder containing full fasta file of current strain
        os.chdir(file_path + "/" + strain)
    
        # open and parse the input FASTA file 
        fasta_input_file = SeqIO.parse(strain + ".fasta", "fasta")
    
        # dict-type variable to store and return accession no. (key)-fasta_record (value) pairs
        fasta_dict = {}
        
        # iterate through FASTA file and record protein accession no., description, sequence, and length as fasta_record object
        for protein in fasta_input_file:
            temp_fasta_obj = fasta_record(protein)
            fasta_dict[temp_fasta_obj.qacc] = temp_fasta_obj
            
        # store parsed dict for current strain in outer dict
        strains_fasta_dict[strain] = fasta_dict

        # should not need to explicitly close fasta_input_file because filename and not handle was used
        #fasta_input_file.close()
    
    return strains_fasta_dict
     
    
#%% Main program
    
# process all strains for master dict
# first item in tuple is list of names of all subject organisms
# outer dict is strains_dict[strain_name]:process_blast_file_dict
# inner dict is process_blast_file_dict[qacc]:blastp_record objects
(subj_strains_list, strains_dict) = process_strains(subject_organism_filename, query_organism, io_subfolder, blast_file_name_format)

# use POI qacc to search against master strains_dict for specific blastp_record objects
# outer dict is POI_blast_dict[qacc]:strains_dict
# inner dict is strains_dict[strain_name]:blastp_record object for specific strain and qacc
POI_blast_dict = search_blast(strains_dict, query_proteins_filename, io_subfolder)

# convert blastp scores for all strains and protein queries into a Pandas DataFrame of float values
blast_df = dict_to_df(POI_blast_dict, score_type, qcovhsp_cutoff)

# obtain list of all strains analyzed (query and subject organisms)
full_strains_list = subj_strains_list.copy()
full_strains_list.append(query_organism)

# process full fasta file for query and all subject organisms into dict
strains_fasta_dict = fasta_to_dict(full_strains_list, io_file_path)


#%% Generate plots

# heatmap of blastp percent identity scores for each protein query and strain
# this heatmap accounts for percent query coverage (i.e. qcovhsp_cutoff)
#clustermap = seaborn.clustermap(blast_df, cmap = "coolwarm", row_cluster = True, col_cluster=False, cbar_kws={"label":score_type})
clustermap = seaborn.clustermap(blast_df, cmap = "coolwarm", row_cluster = True, col_cluster=False, cbar_kws={"label":score_type + " * qcovhsp / 100"})
clustermap.ax_row_dendrogram.set_visible(False)
plt.title("Comparative blast results")
plt.show()

# individual histograms of identity scores for all protein queries per strain
for column_name, column_data in blast_df.iteritems():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
        ax.spines[axis].set_color("gray")
    plt.title(column_name, fontsize=18.0)
    plt.xlabel(score_type, fontsize=18.0)
    plt.ylabel('number of protein queries', fontsize=18.0)
    plt.hist(column_data, bins=20, alpha=0.6, edgecolor="black", linewidth=1.0)
    plt.xlim([0.0,100.0])
    plt.grid(False)
    plt.tick_params(axis="both", direction="out", length=5.0, width=2.0, color="gray", labelsize=18.0)
    plt.show()

# overlayed histograms of all identity scores per strain
blast_df.plot.hist(bins=20, alpha=0.5, edgecolor="black")
plt.title("All strains overlayed")
plt.xlabel(score_type)
plt.ylabel('number of protein queries')
plt.xlim([0.0,100.0])
plt.grid(False)
plt.show()

#%% Print out files

# # indicate gene for which homologs will be extracted and written to out file
# protein_of_interest = "GGFLHMPP_01381"

# # navigate to path to print out file
# os.chdir(io_subfolder)
# out_file = open(protein_of_interest + "_homologs.fasta", "w")

# # iterate through all subject organisms
# # extract top-scoring homolog (top hit, top HSP) of gene of interest in each subject organism
# # use homolog accession no to retrieve its amino acid sequence
# for strain, blastp_record_object in POI_blast_dict[protein_of_interest].items():
    
#     homolog_saccver = blastp_record_object.hsp_dict["saccver"]
#     out_file.write(">" + strain + "_" + homolog_saccver + "\n")
    
#     fasta_dict = strains_fasta_dict[strain]
#     strain_fasta_rec = fasta_dict[homolog_saccver]
    
#     out_file.write(strain_fasta_rec.seq + "\n")

# # finally print out query organism gene of interest sequence
# out_file.write(">" + query_organism + "_" + protein_of_interest + "\n")
# out_file.write(strains_fasta_dict[query_organism][protein_of_interest].seq)

# out_file.close()

# out_file_2 = open("EcFhuA-BcepHomologs.fasta", "w")
# for gene in POI_blast_dict.keys():
    
#     out_file_2.write(">" + gene + "\n")
#     out_file_2.write(strains_fasta_dict[query_organism][gene].seq + "\n")
# out_file_2.close()
    


    