#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 14:39:35 2020

@author: trucdo
"""

#%% Define path, file names, and other parameters

# indicate name (exclude file extension) and path of fasta-formatted .txt file containing TFBS
# output will save to a genome_file_name_PWMmodel subfolder in this path
in_path = "path_to_directory"
tfbs_file_name = "Paeruginosa_fur_sites"

# indicate name of Genbank file (exclude file extension) of organism to scan for TFBS
genome_file_name = "name_of_Genbank_file_without_extension"

# indicate background frequency of each nucleotide for the organism under investigation
nucleotide_prob_dict = {"A": 0.165, "G": 0.335, "C": 0.335, "T": 0.165}

# indicate percent cutoff of top hits
percent_cutoff = 5

#%% Import modules

# to interact with OS-dependent functionality
import os

# to store transcription factor binding site as a Pandas dataframe and to work with arrays
import pandas as pd
#import numpy as np

# to handle Genbank and FASTA files with SeqIO subpackage from Bio package
from Bio import SeqIO

# for some mathematical calculations
import math

# to graph sequence logo
import logomaker as lm

# import matplotlib package for plotting data
import matplotlib.pyplot as plt

# to calculate and record process time
import time

# to handle regex search
import re

#%% Function definitions

def fasta_to_df(in_path, fasta_file_name):
# function to convert known TFBSs formatted as a fasta file into a pandas dataframe
# inputs: path to fasta file (str); name of fasta file (str)
# output: tuple of number of TFBSs and dataframe summarizing TFBSs
# in the dataframe, each row is a different TFBS and each column is a nucleotide within the site
    
    # navigate to folder containing fasta-formatted file of TFBSs
    os.chdir(in_path)
    
    # store each binding site as a dict of dict
    # key is name of site and value is an inner dict
    binding_site_dict = {}
    
    # open and parse the input FASTA file 
    fasta_input_file = SeqIO.parse(fasta_file_name + ".txt", "fasta")
    
    # step through each binding site
    # user counter to keep track of number of TFBS
    
    no_tfbs = 0

    for binding_site in fasta_input_file:
      
        # inner dict in which key is the nucleotide position in the sequence and value is the nucleotide
        sequence_dict = {}
        
        # retrieve position and identity of the nucleotide at the position (i.e. parse the seq into a dict)
        # turn all nucleotides to uppercase
        for position in range(len(binding_site.seq)):
            sequence_dict[position+1] = binding_site.seq[position].upper()
        
        # assign parsed seq to dict with the key as TFBS name
        binding_site_dict[binding_site.description] = sequence_dict
        
        no_tfbs = no_tfbs + 1

    # close input fasta file handle
    # comment out beecause Della gives an error saying that fasta_input_file has no attribute close
    #fasta_input_file.close()
    
    # convert dict to dataframe and transpose so each binding site is a row and each column is a nucleotide in the seq
    binding_site_dataframe = pd.DataFrame(binding_site_dict).transpose()
    
    return (no_tfbs, binding_site_dataframe)


def df_to_pfm(tfbs_df):
# function to convert pandas dataframe summmarizing TFBSs to a position frequency matrix
# input: dataframe (output from fasta_to_df function)
# output: dataframe that summarizes the number and type of nucleotides at each position across all TFBSs
    
    # dict to store frequency of each nucleotide at each position across all TFBSs
    # key is position in the TFBS and value is inner dict
    frequency_dict = {}

    # iterate through the columns (i.e. position in seq) of the dataframe
    for position in tfbs_df.columns:
        
        # extract all nucleotides across all TFBS at current specific position in loop
        # then concatenate nucleotides in the list into a string
        nucl_list = tfbs_df[position].to_list()
        nucl_str = "".join(nucl_list)
             
        # inner dict wherein key is nucleotide and value is frequency of the nucleotide across all TFBS at a position
        nucleotide_dict = {}
        nucleotide_dict["A"] = nucl_str.count("A")
        nucleotide_dict["C"] = nucl_str.count("C")
        nucleotide_dict["G"] = nucl_str.count("G")
        nucleotide_dict["T"] = nucl_str.count("T")
        nucleotide_dict["total"] = len(nucl_str)
        
        frequency_dict[position] = nucleotide_dict
        
    # convert dict to dataframe
    pfm_df = pd.DataFrame(frequency_dict)
        
    return pfm_df


def pfm_to_pwm(tfbs_pfm, no_tfbs, nucl_prob_dict):
# function to convert position frequency matrix of TFBS to position weight matrix
# inputs: dataframe (index 1 output from df_to_pfm function); number of TFBSs (int; index 0 output of df_to_pfm)
# inputs: dict of background probability for each nucleotide (key is str; value is int)
# output: tuple in which
# index 0: dataframe that summarizes the log-scaled normalized frequency value of each nucleotide at each TFBS position
# index 1: information content (float) of calcluated position weight matrix

    # dict to store log-scaled pseudocount-corrected frequency at each position for all TFBSs
    # key is position in the TFBS and value is inner dict
    weight_dict = {}
    
    # to store information content of PWM
    info_content = 0

    # iterate through the columns (i.e. position in seq) of the dataframe
    for position in tfbs_pfm.columns:
        
        # extract the data series corresponding to the current column
        column_series = tfbs_pfm[position]
        
        # extract the total number of nucleotides at current position (equals no_tfbs in most cases)
        no_nucl = column_series["total"]
        
        # calculate pseudocount as square root of number of TFBSs used to build the model divided by 4 for the 4 nucl
        pseudocount = math.sqrt(no_tfbs) / 4
        
        # inner dict wherein key is nucleotide and value is weighted frequency of the nucleotide at current position
        nucleotide_dict = {}
        
        # step through each nucleotide
        for nucleotide in ["A", "C", "G", "T"]:
            
            # extract the raw count of current nucleotide in the current column
            nucleotide_count = column_series.loc[nucleotide]
            
            # calculate corrected probability for current nucleotide in current column
            cor_prob = (nucleotide_count + pseudocount) / (no_nucl + (4 * pseudocount))
            
            # calculate position weight matrix value of current nucleotide in current position
            # divide corrected nucleotide probability by expected background probability of current nucleotide
            # convert to log base 2
            nucleotide_dict[nucleotide] = math.log(cor_prob/nucl_prob_dict[nucleotide], 2)
            
            # update information content
            info_content = cor_prob * nucleotide_dict[nucleotide]
            
        weight_dict[position] = nucleotide_dict
        
    # convert dict to dataframe
    pwm_df = pd.DataFrame(weight_dict)

    return (pwm_df, info_content)


def find_consensus_seq(tfbs_pwm):
# function to define the highest-scoring sequence (i.e. consensus sequence) and max score from a position weight matrix
# input: position weight matrix dataframe (index 0 output of pfm_to_pwm function)
# output: tuple of consensus sequence (index 0) and its score (index 1)

    # initiate variables to store consensus sequence (str) and score (float)
    consensus_score = 0
    consensus_seq = ""

    # iterate through the columns (i.e. position in seq) of the dataframe
    for position in tfbs_pwm.columns:
        
        # extract the highest nucleotide PWM score at current position
        consensus_score = consensus_score + max(tfbs_pwm[position])
        
        # create a bool filter to identity best nucleotide at current position
        max_filter = tfbs_pwm[position] == max(tfbs_pwm[position])
        
        # apply bool filter to current column to extract identity of highest-scoring nucleotide
        # index method produces pandas index object, so convert to list object
        max_nucl_list = list(tfbs_pwm[position][max_filter].index)
        
        # in case there are multiple nucleotides that are equally high-scoring at current position
        # store both equally-scoring nucleotides separated by a pipe character
        # and put parentheses around the possible alternate nucleotides at this position
        # if max_nucl_list has only 1 item (i.e. only one high-scoring nucleotide), then no pipe character is added
        if len(max_nucl_list ) > 1:
            max_nucl = "(" + "|".join(max_nucl_list) + ")"
        else:
            max_nucl = "".join(max_nucl_list)
        
        consensus_seq = consensus_seq + max_nucl
        
    return (consensus_seq, consensus_score)
        

def seq_score(tfbs_pwm, dna_seq):
# function to score a DNA sequence based on its degree of match to a TFBS described by a position weight matrix
# inputs: position weight matrix dataframe (output of pfm_to_pwm function); DNA sequence of interest (str)
# output: dataframe summarizing weighted score for each nucleotide of the DNA sequence
# each column is a nucleotide and there are two rows
# the last column is the total score (sum of all nucleotide PWM scores), complete DNA and TFBS sequence
# the first row is the scores and the second row is the nucleotide/DNAsequence

    # check that the DNA sequence has the same length as the TFBSs and return error if different lengths
    
    if len(dna_seq) != len(tfbs_pwm.columns):
        print("seq_score function failed: DNA sequence length differs from TFBS length")
        return "seq_score function failed: DNA sequence length differs from TFBS length"
    else:
        
        # ensure DNA query sequence is all uppercase
        DNA_seq = dna_seq.upper()
        
        # dict to store PWM score of each nucleotide in the DNA sequence for all positions
        # key is position in the TFBS and value is inner dict of nucleotide scores/identity
        score_dict = {} 
        
        # iterate through each nucleotide of the DNA sequence
        # keep track of the nucleotide position to query the correct PWM column
        column_no = 1
        # keep track of sum score of all nucleotides
        sum_score = 0;
        
        for nucleotide in DNA_seq:
            
            # inner dict wherein key is nucleotide "score" and "identity" and values follow
            inner_score_dict = {}
            
            # extract PWM value
            inner_score_dict["score"] = tfbs_pwm.loc[nucleotide,column_no]
            
            # extract current nucleotide
            inner_score_dict["identity"] = nucleotide
                                
            # store inner dict in outer score_dict
            score_dict[column_no] = inner_score_dict
            
            # update sum_score
            sum_score = sum_score + inner_score_dict["score"]
            
            # increase column count by 1 to shift to next column
            column_no = column_no + 1
            
        # convert score dict to dataframe
        score_df = pd.DataFrame(score_dict, index = ["score", "identity"])
        
        # add a new final column to store sum of all nucleotide scores and complete DNA sequence
        score_df["sum"] = [sum_score, DNA_seq]
        
        return score_df
    

def gbk_to_fasta(in_path, genome_file_name):
# function to convert whole-genome Genbank seqence file for query organism to a FASTA file
# calculates both the forward and reverse-complement of genome sequence
# inputs: path to genome file (str); name of genome file (str)
# output: returns dict of seq_record objects; new fasta file is saved to in_path
# dict keys are names of chromosomes (e.g. chr_1)
# dict values are the seq_record object corresponding to the chromosome

    # navigate to folder containing genome sequence file
    os.chdir(in_path)
    
    # open Genbank file to read 
    input_handle  = open(genome_file_name + ".gbk", "r")
    output_handle = open(genome_file_name + ".fasta", "w")
    
    # extract DNA sequence of all chromosomes and store in FASTA format
    # for each chromosome, also write out the reverse complement sequence
    # each FASTA entry is a different chromosome
    # store each seq_record (forward chromosome) as a key-value pair in dict
    
    seq_record_dict = {}
    
    for seq_record in SeqIO.parse(input_handle, "genbank"):
        
        # store seq_record object in dict
        seq_record_dict["chr_" + str(seq_record.name)] = seq_record
    
        # sequence as in Genbank record
        forward_seq = str(seq_record.seq)
        
        # reverse-complement of sequence
        reverse_seq = forward_seq.replace("A", "t")
        reverse_seq = reverse_seq.replace("C", "g")
        reverse_seq = reverse_seq.replace("G", "c")
        reverse_seq = reverse_seq.replace("T", "a")
        reverse_seq = reverse_seq.upper()
        reverse_seq = reverse_seq[::-1]
        
        output_handle.write(">chr_" + str(seq_record.name) + "\n" + forward_seq+ "\n")
        output_handle.write(">chr_" + str(seq_record.name) + "_reverse-complement\n" + reverse_seq+ "\n")

    output_handle.close()
    input_handle.close()
    
    return seq_record_dict    

    
def scan_genome(in_path, genome_fasta_file, tfbs_pwm):
# function to scan whole-genome FASTA sequence file for candidate TFBSs using position weight matrix
# inputs: path to genome file (str); name of FASTA-formatted genome file (str; file output of gbk_to_fasta)
# inputs: position weight matrix dataframe (output of pfm_to_pwm function)
# output: dataframe summarizing sequence and genome position of all sequences with length equal to TFBS length
# dataframe also contains match score calculated for each sequence using PWM

    # navigate to folder containing genome sequence file
    os.chdir(in_path)
    
    # open FASTA file to read 
    input_handle  = open(genome_fasta_file + ".fasta", "r")
      
    # length of TFBS
    tfbs_len = len(tfbs_pwm.columns)
    
    # outer dict to store extracted sequence
    # key is arbitrary number keeping track of order when sequence was extracted
    scan_dict = {}
    count = 1
    
    # scan each chromosome and extract sequence of the same length as TFBS
    for chromosome in SeqIO.parse(input_handle, "fasta"):
        
        # tack on the first X nucleotides of the chromosome to the end of the chromosome to account for circular chromosome
        # where X is the length of TFBS minus 1   
        cir_chr = str(chromosome.seq).rstrip("\n") + str(chromosome.seq[0:(tfbs_len-1)]).rstrip("\n")
        cir_chr = cir_chr.upper()
        
        # use a sliding window of 1 base and extract sequence and its location on the chromosome
        # scan to the end of the chromosome
        for i in range(0, len(cir_chr) - tfbs_len + 1):
            
            # inner dict wherein keys: sequence, chromosome_no, location
            # all values are strings except for "location" key, whose value is an int
            inner_dict = {}
            
            inner_dict["sequence"] = cir_chr[i:i+tfbs_len]
            inner_dict["chromosome_no"] = str(chromosome.name)
            inner_dict["start_location"] = i+1
            
            # calculate match between the extracted sequence and TFBS using PWM
            # seq_score output is a dataframe of nucleotide score and identity
            # just extract the total PWM score for entire query sequence
            inner_dict["PWM_score"] = seq_score(tfbs_pwm, inner_dict["sequence"]).loc["score","sum"]
                    
            # store inner_dict for current sequence
            scan_dict[count] = inner_dict
            
            # update count before proceeding to next sequence
            count = count + 1
        
    # convert dict to dataframe and transpose so each extracted sequence is a row and its info is a column
    extracted_fragments = pd.DataFrame(scan_dict).transpose()
    
    return extracted_fragments


def filter_cds(seq_record_dict, fragments_df):
# function to filter TFBS matches for those entirely within non-coding (i.e. regulatory) regions
# inputs: dict of SeqRecord objects (return output of gbk_to_fasta function)
# inputs: dataframe of fragments to assess (output of scan_genome before or after PWM filtering)
# output: filtered dataframe of all sequence fragments that are fully wthin non-coding regions
# dataframe also contains genome position and match score calculated for each sequence using PWM

    # iterate through all key-value pairs (i.e. all chromosomes and their SeqRecord objects)
    # add range of coordinates of all CDS to a set stored in a dict under chromosome name key
    # note that all coordinates correspond to forward-strand position
    cds_coor_dict = {}
    for chr_name, seq_record in seq_record_dict.items():
        # set to store coordinates of all CDS for current chromosome
        cds_coor_set = set()
        for seq_feature in seq_record.features:
            if seq_feature.type=="CDS":
                # account for the start and end location numbering following Python's zero-index scheme
                fstart = int(seq_feature.location.start) + 1
                fend = int(seq_feature.location.end)
                # add coordinates within the CDS range to set of all CDS coordinates in current chromosome
                cds_coor_set.update(list(range(fstart, fend + 1)))
        # store set of CDS coordinates for current chromosome in dict        
        cds_coor_dict[chr_name] = cds_coor_set

    # iterate through all the sequence fragments (each row) stored in dataframe
    # check if the entire fragment is present within non-coding regions of the genome
    # add fully non-coding fragments to filtered dict
    
    noncds_dict = {}
    
    for index, row in fragments_df.iterrows():
            
        # retrieve name of current chromosome and set rc_flag to default False for forward strand
        chr_no = row["chromosome_no"]
        rc_flag = False
        
        # if current fragment is found on reverse strand
        # remove the reverse complement from chromosome name for later search
        # flag to convert reverse-complement coordinates to forward coordinates
        if re.search("_reverse-complement", chr_no):
            chr_no = chr_no.rstrip("_reverse-complement")
            rc_flag = True
        
        # retrieve length of current chromosome
        chr_len = len(seq_record_dict[chr_no].seq)
        
        # for forward strand
        # retrieve coordinates as is of all nucleotides that comprise the current sequence fragment as a list
        if rc_flag == False:
            start = row["start_location"]
            interval = len(row["sequence"])
            all_coor_list = list(range(start, start + interval))
        # for reverse-complement strand
        else:
            start = chr_len - row["start_location"] + 1
            interval = len(row["sequence"])
            all_coor_list = list(range(start, start - interval, -1))
        
        # iterate through all coordinates in list and check if nucleotide position is within CDS set for current chr
        # default is fragment fully within non-coding region
        # break search to exit loop if any nucleotide is within CDS set
        cds_flag = False
        
        # to account for extra bases at the end due to wrapping around to position 1 for circular chromosome
        # first convert nucleotides at the end to actual coordinates
        # the last X positions of a chromosome is equal to the first X positions of the same chr
        # where X = TFBS length - 1
        
        for coor in all_coor_list:
            if coor > chr_len:
                coor = coor - chr_len
            if coor in cds_coor_dict[chr_no]:
                cds_flag = True
                break
           
        # if entire sequence fragment is outside CDS, add to dict
        if cds_flag == False:
            noncds_dict[index] = row
            
    # convert dict to dataframe and transpose so each filtered sequence is a row and its info is a column
    noncds_fragments = pd.DataFrame(noncds_dict).transpose()
    
    return noncds_fragments


def find_cds_neighbor(seq_record_dict, fragments_df):
# function to find the adjacent CDS to the candidate TFBS
# inputs: dict of SeqRecord objects (return output of gbk_to_fasta function)
# inputs: dataframe of fragments to assess (output of scan_genome before or after filtering)
# output: dataframe of input sequence fragments with new column describing adjacent CDS

    # iterate through all key-value pairs (i.e. all chromosomes and their SeqRecord objects)
    # note that all coordinates correspond to forward-strand position
    
    # outer dict key is chromosome name and value is multi-layered dict
    chr_dict = {}

    for chr_name, seq_record in seq_record_dict.items():
         
        # middle dict separating start/end coordinates of forward and reverse-complement strands
        strand_dict = {}
        
        # in-most dicts wherein keys are start/end coordinates and values are tuple (CDS locus tags, annotated product)
        # separate start/end and forward/reverse strands into different dicts
        scoor_dict = {}
        ecoor_dict = {}
        rscoor_dict = {}
        recoor_dict = {}
            
        for seq_feature in seq_record.features:
            
            if seq_feature.type=="CDS":     
                
                # extract start/end positions of CDS
                # account for Python vs. Genbank indexing by adding 1 to start position
                start_pos = int(seq_feature.location.start) + 1
                end_pos = int(seq_feature.location.end)
                
                # for forward CDS sequence
                if seq_feature.location.strand == 1:
                    scoor_dict[start_pos] = (seq_feature.qualifiers["locus_tag"][0], seq_feature.qualifiers["product"][0])
                    strand_dict["forward_start"] = scoor_dict
                    ecoor_dict[end_pos] = (seq_feature.qualifiers["locus_tag"][0], seq_feature.qualifiers["product"][0])
                    strand_dict["forward_end"] = ecoor_dict
                # for reverse-complement CDS sequence
                else:
                    rscoor_dict[start_pos] = (seq_feature.qualifiers["locus_tag"][0], seq_feature.qualifiers["product"][0])
                    strand_dict["reverse_start"] = rscoor_dict
                    recoor_dict[end_pos] = (seq_feature.qualifiers["locus_tag"][0], seq_feature.qualifiers["product"][0])
                    strand_dict["reverse_end"] = recoor_dict
        
        # store coordinates of current chromosome in multi-layered dict
        chr_dict[chr_name] = strand_dict
        
            
    # iterate through all candidate TFBS fragments and search for closest CDS on same strand
    cds_info_dict = {}
    
    for index, row in fragments_df.iterrows():
        
        # retrieve name and start position of current chromosome
        # set rc_flag to default False for forward strand
        chr_no = row["chromosome_no"]
        rc_flag = False
        fragment_loc = row["start_location"]
                
        # if current fragment is found on reverse strand
        # remove the reverse complement from chromosome name for later search
        # flag to convert reverse-complement coordinates to forward coordinates
        # convert start position on reverse-complement strand to corresponding forward position strand
        if re.search("_reverse-complement", chr_no):
            chr_no = chr_no.rstrip("_reverse-complement")
            rc_flag = True
            fragment_loc = len(seq_record_dict[chr_no].seq) - row["start_location"] + 1
        
        # to search on forward strand
        if rc_flag == False:
            strand = chr_dict[chr_no]
            coor = strand["forward_start"]
            all_coor = list(coor.keys())
            # add start position of current TFBS fragment to list and change to ascending order
            all_coor.append(fragment_loc)
            all_coor.sort()
            # find index of TFBS fragment in list and retrieve number in index + 1 position (i.e. closest forward CDS)
            closest_CDS_index = all_coor.index(fragment_loc) + 1
            # if TFBS is at the end of the list, then loop back to the first item in list
            if closest_CDS_index > len(all_coor) - 1:
                closest_CDS_index = 0
            # use retrieved position to obtain CDS locus tag and annnotation from multi-layered dict
            closest_CDS_pos = all_coor[closest_CDS_index]
            cds_name = coor[closest_CDS_pos][0]
            cds_annotation = coor[closest_CDS_pos][1]
        # to search on reverse strand
        else:
            strand = chr_dict[chr_no]
            coor = strand["reverse_end"]
            all_coor = list(coor.keys())
            # add start position of current TFBS fragment to list and change to ascending order
            all_coor.append(fragment_loc)
            all_coor.sort()
            # find index of TFBS fragment in list and retrieve number in index - 1 position (i.e. closest reverse CDS)
            # if TFBS is at the beginnning of the list, then loop forward to the last item in list (index of -1)
            closest_CDS_index = all_coor.index(fragment_loc) - 1
            # use retrieved position to obtain CDS locus tag from multi-layered dict
            closest_CDS_pos = all_coor[closest_CDS_index]
            cds_name = coor[closest_CDS_pos][0]
            cds_annotation = coor[closest_CDS_pos][1]
            
        # convert current row (i.e. current sequence fragment info) to dict
        # add two new key-value pairs for CDS_neighbor and CDS_annotation
        temp_dict = dict(row)
        temp_dict["CDS_neighbor"] = cds_name
        temp_dict["CDS_annotation"] = cds_annotation
        cds_info_dict[index] = temp_dict
        
    # convert dict to dataframe and transpose so each filtered sequence is a row and its info is a column
    cds_info = pd.DataFrame(cds_info_dict).transpose()
    
    return cds_info
        
    
#%% Main program

# record and print the process start time in second
start_time = time.time()
print("Start time: " + str(start_time) + " seconds")

# convert TFBSs into dataframe
no_tfbs, binding_site_df = fasta_to_df(in_path, tfbs_file_name)
# convert TFBS dataframe into position frequency matrix
tfbs_pfm = df_to_pfm(binding_site_df)
# convert position frequency matrix to position-specific scoring matrix and information content of PWM
tfbs_pwm, info_content = pfm_to_pwm(tfbs_pfm, no_tfbs, nucleotide_prob_dict)
# determine consensus sequence and max score based on position frequency matrix of input TFBSs
consensus_seq, consensus_score = find_consensus_seq(tfbs_pwm)

# convert whole-genome Genbank sequence file of query organism to FASTA format and save dict of SeqRecord objects
seq_record_dict = gbk_to_fasta(in_path,genome_file_name)

# scan entire genome and extract sequence fragments with length equal to TFBS length
genome_scan_df = scan_genome(in_path, genome_file_name, tfbs_pwm)

# filter genome for only those sequences with PWM scores > information content
ICfiltered_genome_df = genome_scan_df[genome_scan_df["PWM_score"] > info_content]

# filter IC-filtered genome further for TFBS matches within non-coding regions
CDS_IC_filtered_genome_df = filter_cds(seq_record_dict, ICfiltered_genome_df)

# also filter entire genome for TFBS matches within non-coding regions
CDS_filtered_genome_df = filter_cds(seq_record_dict, genome_scan_df)

# annotate CDS-IC filtered TFBS matches to find closest CDS
CDS_IC_neighbor_df = find_cds_neighbor(seq_record_dict, CDS_IC_filtered_genome_df)

# also annotate CDS filtered TFBS matches to find closest CDS
CDS_neighbor_df = find_cds_neighbor(seq_record_dict, CDS_filtered_genome_df)

# filter genome for only those sequences with PWM scores > user-defined cutoff
index_cutoff = math.floor(genome_scan_df.shape[0] * (100-percent_cutoff)/100)
ordered_df = genome_scan_df.sort_values("PWM_score", axis = 0, ascending = True)
ufiltered_genome_df = ordered_df.iloc[index_cutoff:,:]

# record and print the process end time in seconds
end_time = time.time()
print("End time: " + str(end_time) + " seconds")

# record and print the process time in min and second
elapsed_time = end_time-start_time
# if the process took less than a minute
if elapsed_time < 60:
    elapsed_min = 0
    elapsed_sec = elapsed_time
# if the process took more than a minute
else:
    elapsed_min = math.floor(elapsed_time/60)
    elapsed_sec = (elapsed_time/60 - elapsed_min) * 60

print("Elapsed time: " + str(elapsed_min) + " minutes and " + str(elapsed_sec) + " seconds")


#%% Data visualization

# make and navigate to new directory to store results
out_path = in_path + "/" + genome_file_name + "_PWMmodel/"
if os.path.exists(out_path):
    os.chdir(out_path)
else:
    os.mkdir(out_path)
    os.chdir(out_path)
        
# plot the sequence logo using the position frequency matrix of input TFBSs
# also save sequence logo to out_path
logo = lm.Logo(tfbs_pwm.transpose(), font_name = 'Arial Rounded MT Bold', fade_below=0.5, shade_below=0.5)
logo.ax.set_xlabel('Position',fontsize=14)
logo.ax.set_ylabel("Position-specific score", labelpad=-1,fontsize=14)
plt.savefig(tfbs_file_name + "_sequencelogo.png", bbox_inches = "tight")

# plot histogram to visualize distribution of PWM scores for all sequences
# save histogram to out_path
fig = plt.figure()
ax = fig.gca()
plt.title("All PWM scores from genome scan", fontsize=18.0)
plt.xlabel("PWM score", fontsize=18.0)
plt.ylabel("Number of sequences", fontsize=18.0)
plt.hist(genome_scan_df["PWM_score"], bins=20, alpha=0.6, edgecolor="black", linewidth=1.0)
# plot vertical line indicating threshold = information content and user-defined cutoff
plt.axvline(x = info_content, ymin = 0, color='seagreen', linewidth = 2)
plt.axvline(x = ufiltered_genome_df.iloc[0]["PWM_score"], ymin = 0, color='orange', linewidth = 2)
plt.xlim([-consensus_score,consensus_score])
plt.grid(False)
plt.tick_params(axis="both", direction="out", length=5.0, width=2.5, color="gray", labelsize=18.0)
ax.spines["bottom"].set_linewidth(1.5)
ax.spines["top"].set_linewidth(1.5)
ax.spines["left"].set_linewidth(1.5)
ax.spines["right"].set_linewidth(1.5)
plt.savefig(genome_file_name + "_predictedtfbs.png", bbox_inches = "tight")
plt.show()


#%% Write out results

tfbs_pfm.to_csv(tfbs_file_name + "_pfm.csv")
tfbs_pwm.to_csv(tfbs_file_name + "_pwm.csv")
genome_scan_df.to_csv(genome_file_name + "_allfragments.csv")
ICfiltered_genome_df.to_csv(genome_file_name + "_ICfilteredfragments.csv")
ufiltered_genome_df.to_csv(genome_file_name + "_Ufilteredfragments.csv")
CDS_IC_neighbor_df.to_csv(genome_file_name + "_CDS-IC-filteredfragments.csv")
CDS_neighbor_df.to_csv(genome_file_name + "_CDSfilteredfragments.csv")

out_file = open(genome_file_name + "_stats.txt", "w")
out_file.write("----------PWMmodel.py----------\n\n")
out_file.write("Input file(s):" + tfbs_file_name + ", " + genome_file_name + "\n")
out_file.write("Output file(s): " + tfbs_file_name + "_pfm.csv, " + tfbs_file_name + "_pwm.csv, ")
out_file.write(genome_file_name + "_allfragments.csv, " + genome_file_name + "_ICfilteredfragments.csv, ")
out_file.write(genome_file_name + "_Ufilteredfragments.csv, " + genome_file_name + "_CDS-IC-filteredfragments.csv, ")
out_file.write(genome_file_name + "_CDSfilteredfragments.csv, " + genome_file_name + "_stats.txt, " + genome_file_name + ".fasta, ")
out_file.write(tfbs_file_name + "_sequencelogo.png, " + genome_file_name + "_predictedtfbs.png\n")
out_file.write("Number of TFBSs used to build model: " + str(no_tfbs) + "\n")
out_file.write("Length of TFBSs: " + str(binding_site_df.shape[1]) + " bases\n")
out_file.write("Background probability of nucleotides: " + str(nucleotide_prob_dict) + "\n")
out_file.write("TFBS Consensus sequence: " + consensus_seq + "\n")
out_file.write("PWM score for perfect match to consensus sequence: " + str(consensus_score) + "\n")
out_file.write("Information content (IC): " + str(info_content) + "\n")
out_file.write("User-defined percent cutoff for top hits: " + str(percent_cutoff) + "\n")
out_file.write("Number of genome fragments assessed: " + str(genome_scan_df.shape[0]) + "\n")
out_file.write("Number of IC-filtered genome fragments: " + str(ICfiltered_genome_df.shape[0]) + "\n")
out_file.write("Number of CDS and IC-filtered genome fragments: " + str(CDS_IC_filtered_genome_df.shape[0]) + "\n")
out_file.write("Number of CDS-filtered genome fragments: " + str(CDS_filtered_genome_df.shape[0]) + "\n")
out_file.write("Number of user-filtered genome fragments " + str(ufiltered_genome_df.shape[0]) + "\n")
out_file.write("Elapsed time: " + str(elapsed_min) + " minutes and " + str(elapsed_sec) + " seconds\n")
out_file.close()
