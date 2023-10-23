#!/usr/bin/env python3

import os
import csv
import pandas
pandas.set_option('display.max_colwidth', None)
from wolfe import *

path = '/mnt/scratch/marshall/comparative_genomics/bacteroides/hmmer/results'
path_to_orfs = '/mnt/scratch/marshall/comparative_genomics/bacteroides/orfs/extracted_orfs'
path_to_genomes = '/mnt/scratch/marshall/comparative_genomics/bacteroides/database/genomes'
output = "all_captured_data.txt"
with_seqs_output = "all_hits_evalue_seqs_testing.csv"

def capture_tsv_info(file_name):
    #this captures all of entries in the tsv output files. It excludes the
    #header which describes each column and the end of the file which contains
    #information about the tsv file
    open_file = pandas.read_csv(file_name, sep = '\t')
    orf_info = open_file.iloc[2:-10,]
    return str(orf_info)

def pull_name(file_name):
    #this pulls the full file name, then splits on / to just catch
    #the final part of the name. It splits again on .faa and takes
    #the first indexed string which I expect to be the species name
    full_path = os.fsdecode(file_name)
    full_file_name = full_path.split('/')[-1]
    species_name = (full_file_name.split('.faa')[0])
    return species_name

def pull_orf(data_frame):
    #using wolfe.py tools, ingests fasta file of ORFS, finds the entry
    #with the matching header, and outputs the sequence
    name_of_species = data_frame['species_name']
    name_of_entry = (data_frame['ORF_name'] + ' ' + data_frame['gene_location']
    + ' ' + data_frame['gene_type'] + ' ' + data_frame['length'] + ' ' +
    data_frame['frame'] + ' ' + data_frame['start'] + ' ' + data_frame['stop'])
    #the name is so verbose so that it exactly matches the fasta header
    infasta = FastaFile()
    with open (path_to_orfs + "/" + name_of_species + '.faa') as inf:
        infasta.read_whole_file(inf)
        out_protein_seq = FastaEntry()
        out_protein_seq  = infasta.pull_entry(name_of_entry)
        return out_protein_seq.pull_seq(0, -1)

def pull_ops(data_frame):
    name_of_species = data_frame['species_name'] 
    strand = data_frame['gene_location'][-2]
    coordinates = data_frame['gene_location'].replace('[', "")
    coordinates = coordinates.replace("]","")
    coordinates = coordinates[0:-3]
    chromosome = data_frame['ORF_name'].split("_ORF.")
    chromosome = chromosome[0]
    if strand == '+':
        start_plus = int(coordinates.split('-')[0])
        end_plus = start_plus-200
        infasta = FastaFile()
        with open (path_to_genomes + "/" + name_of_species) as inf:
            infasta.read_whole_file(inf)
            for header in infasta.chrm_names():
                try:
                    if header.split(" ")[0] == chromosome:
                        full_seq = infasta.pull_entry(header)
                        ops_seq_plus = full_seq.pull_seq(end_plus,start_plus)
                        return ops_seq_plus
                except ValueError:
                    print("%s, coordinates %s must be circular for result, plus strand"%(name_of_species,coordinates))
                    return 'Too_close_to_end'
                    continue

    elif strand == '-':
        start_minus = int(coordinates.split('-')[1])
        end_minus = start_minus+200
        infasta = FastaFile()
        with open (path_to_genomes + "/" + name_of_species) as inf:
            infasta.read_whole_file(inf)
            for header in infasta.chrm_names():
                try:
                    if header.split(" ")[0] == chromosome:
                        full_seq = infasta.pull_entry(header)
                        ops_seq_plus = full_seq.pull_seq(start_minus,end_minus)
                        return ops_seq_plus
                except ValueError:
                    print("%s, coordinates %s must be circular for result, minus strand"%(name_of_species,coordinates))
                    return 'Too_close_to_end'
                    continue




for result in os.scandir(path):
    name = pull_name(result) #pull_name catches the species name
    hmmer_output = capture_tsv_info(result) 
    each_line = hmmer_output.splitlines() #splits each line of the hmmer_output
                                          #into its own list
                                          #for each orf entry, pulls info
    try:
        for line in each_line[1:]:            
            split_line = line.split()
            target_name = split_line[1]
            e_value = split_line[5]
            gene_location = split_line[19]
            gene_type = split_line[20]
            length = split_line[21]
            frame = split_line[22]
            start = split_line[23]
            stop = split_line[24]

            f = open(output, "a")
            f.write(name + ',' + target_name + ',' + e_value + ',' + gene_location
                    + ',' + gene_type + ',' + length + ',' + frame + ',' + start + 
                    ',' + stop + '\n' )
            f.close
    except IndexError:          
        print("%s had zero identified proteins with homology"%(result))
        continue 
        #Species that have no homologs identified will throw an index error.
        #The error gets handled here and those species with no hits never get
        #captured in any files. Counting the number of things that get thrown
        #out with this error will be the number of species that had no hits
        
captured_hmmer_output = pandas.read_csv(output)
column = ['species_name','ORF_name','e_value','gene_location','gene_type','length','frame','start','stop']
captured_hmmer_output.columns = column

protein_seq = (captured_hmmer_output.apply(pull_orf, axis=1))
captured_hmmer_output['protein_seq'] = protein_seq
ops_seq = (captured_hmmer_output.apply(pull_ops, axis=1))
captured_hmmer_output['ops_seq'] = ops_seq

captured_hmmer_output.to_csv(with_seqs_output)

os.remove(output)

