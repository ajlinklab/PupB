#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 16:22:22 2020

@author: trucdo
"""

#%% Define path, file names, and other parameters

# indicate name (exclude file extension) and path of .csv file containing TFBS hits
in_path = "path_to_directory"
tfbs_file_name = "Burkholderia_cepacia_ATCC_25416_CDS-IC-filteredfragments"

# indicate information content or other value for cutoff (float)
cutoff = 0

#%% Import modules

# to interact with OS-dependent functionality
import os

# to store transcription factor binding site as a Pandas dataframe and to work with arrays
import pandas as pd
#import numpy as np

# import matplotlib and seaborn packages for plotting data
import matplotlib.pyplot as plt

#%% Class definitions

class tfbs(object):
# class to store information for a TFBS
# input: a single TFBS hit from PWMmodel.py

    # constructor to call an instance of the class and parse input to calculate attributes
    # input: self variable referring to PWMmodel.py result for a single TFBS hit (str)
    # all attributes are string type
    
    def __init__(self, full_result):
        
        # split full TFBS result into a list using comma delimiter ","
        parsed_result = full_result.split(",")
        
        self.index = parsed_result[0]
        self.sequence = parsed_result[1]
        self.chromosome_no = parsed_result[2]
        self.start_location = parsed_result[3]
        self.PWM_score = parsed_result[4]
        self.CDS_neighbor = parsed_result[5]
        self.CDS_annotation = parsed_result[6]


#%% Main program

# navigate to directory where results are stored
os.chdir(in_path)

# open TFBS results .csv file, read in each TFBS as a separate line, remove header
tfbs_file_handle = open(tfbs_file_name + ".csv")
tfbs_file_readlines = tfbs_file_handle.readlines()
header = tfbs_file_readlines.pop(0)

# iterate through each TFBS, convert it to a tfbs object, and then dict of dict
# dict key is arbitrary index 
# dict values are the original index, sequence, chromosome_no, start_location, PWM_score of TFBS hit, CDS info
# use inner dict to separately store attributes of each TFBS

all_tfbs_dict = {}
count = 1

out_file = open(tfbs_file_name + "_scores.txt", "w")

for line in tfbs_file_readlines:
    
    inner_dict = {}
    
    tfbs_obj = tfbs(line.rstrip("\n"))
    
    inner_dict["index"] = tfbs_obj.index
    inner_dict["sequence"] = tfbs_obj.sequence
    inner_dict["chromosome_no"] = tfbs_obj.chromosome_no
    inner_dict["start_location"] = int(tfbs_obj.start_location)
    inner_dict["PWM_score"] = float(tfbs_obj.PWM_score)
    inner_dict["CDS_neighbor"] = tfbs_obj.CDS_neighbor
    inner_dict["CDS_annotation"] = tfbs_obj.CDS_annotation
    
    all_tfbs_dict[count] = inner_dict
    
    count = count + 1
    
    out_file.write(str(tfbs_obj.PWM_score + "\n"))
    
out_file.close()
    
# convert TFBS dict to dataframe and transpose so each extracted sequence is a row and its info is a column
all_tfbs_df = pd.DataFrame(all_tfbs_dict).transpose()

# filter TFBS dataframe for only those TFBS hits with PWM > IC cutoff
filtered_tfbs_df = all_tfbs_df[all_tfbs_df["PWM_score"] > cutoff]

#%% Data visualization

# plot histogram to visualize distribution of all PWM scores
# save histogram to in_path
fig = plt.figure()
ax = fig.gca()
plt.title("PWM scores from all fragments", fontsize=18.0)
plt.xlabel("PWM score", fontsize=18.0)
plt.ylabel("Number of sequences", fontsize=18.0)
plt.hist(all_tfbs_df["PWM_score"], bins=20, alpha=0.6, edgecolor="black", linewidth=1.0)
plt.xlim([min(all_tfbs_df["PWM_score"]),max(all_tfbs_df["PWM_score"])])
plt.grid(False)
plt.tick_params(axis="both", direction="out", length=5.0, width=2.5, color="black", labelsize=18.0)
ax.spines["bottom"].set_linewidth(1.5)
ax.spines["top"].set_linewidth(1.5)
ax.spines["left"].set_linewidth(1.5)
ax.spines["right"].set_linewidth(1.5)
plt.savefig(tfbs_file_name  + "_allfragments_histogram.png", bbox_inches = "tight")
plt.show()

# plot histogram to visualize distribution of PWM scores for filtered sequences
# save histogram to in_path
fig = plt.figure()
ax = fig.gca()
plt.title("PWM scores from cutoff-filtered hits", fontsize=18.0)
plt.xlabel("PWM score", fontsize=18.0)
plt.ylabel("Number of sequences", fontsize=18.0)
plt.hist(filtered_tfbs_df["PWM_score"], bins=20, alpha=0.6, edgecolor="black", linewidth=1.0)
plt.xlim([cutoff,max(filtered_tfbs_df["PWM_score"])])
plt.grid(False)
plt.tick_params(axis="both", direction="out", length=5.0, width=2.5, color="black", labelsize=18.0)
ax.spines["bottom"].set_linewidth(1.5)
ax.spines["top"].set_linewidth(1.5)
ax.spines["left"].set_linewidth(1.5)
ax.spines["right"].set_linewidth(1.5)
plt.savefig(tfbs_file_name  + "_cutoff_histogram.png", bbox_inches = "tight")
plt.show()    