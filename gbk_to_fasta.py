#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 18:13:30 2020

@author: trucdo

"""
# this script converts a Genbank file to FASTA format with description recorded but without chromosome number
# source code: https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank2fasta/
# Define: current working directory; gbk_filename; fasta_filename

# allows user to interact with OS-dependent functionality 
import os

# change the current working directory
os.chdir("path_to_directory")

# import SeqIO subpackage from Bio package
from Bio import SeqIO

# indicate name of Genbank file
gbk_filename = "name_of_Genbank_file_with_extension"

# indicate name of FASTA file
fasta_filename = "name_of_FASTA_file_with_extension"

# open Genbank file to read and FASTA file to write
input_handle  = open(gbk_filename, "r")
output_handle = open(fasta_filename, "w")

# step through a SeqRecord iterator
for seq_record in SeqIO.parse(input_handle, "genbank") :
    print ("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            # I added line 40 to the original source code to handle Genbank files containing CDS with missing translation...
            # ...which throws a KeyError because there is no "translation" key in the qualifiers dict within the seq_feature...
            # ...object. So CDSs with missing translation are excluded from the final FASTA file.
            if seq_feature.qualifiers.get("translation", "no translation") != "no translation":
                assert len(seq_feature.qualifiers['translation'])==1
                output_handle.write(">%s %s\n%s\n" % (
                       seq_feature.qualifiers['locus_tag'][0],
                       seq_feature.qualifiers['product'][0],
                       seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()