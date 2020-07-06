#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""closest_genes.py
    Usage: closest_genes.py -f <FILE> -a <FILE> -t <FILE> -m <mode> [-k <INT>]

    Options:
        -h, --help                  show this
        -f, --file <FILE>           gff or bed file with gene coordinates
        -a, --ai_file <FILE>        Alienness index full file
        -t, --tree_results <FILE>   Tree results file
        -m, --mode <mode>           mode to run
        -k, --k_value <INT>         number of close genes to process [default: 5]

    Arguments:
        mode: Eugene   - gff from Eugene Pipeline
              Augustus - gff from Augusts/BRAKER Pipeline
              Bed      - user created bed file
"""

import os
import sys
import re
import math
from docopt import docopt
from pybedtools import BedTool

def main():
    args = docopt(__doc__)
    genome_annotation_file = args['--file']
    ai_file = args['--ai_file']
    tree_results_file = args['--tree_results']
    k_value = int(args['--k_value'])
    mode = args['--mode']

    if mode == "Bed":
        ann_file = BedTool(genome_annotation_file)
    elif mode == "Eugene" or mode == "Augustus":
        ann_file = BedTool(gff_2_bed(genome_annotation_file,mode),from_string=True)
    else:
        sys.exit("\t[!] Specify a valid mode (Bed | Eugene | Augustus)")

    ai_file_in = open(ai_file, 'r')
    gene_classification = {}
    for line in ai_file_in:
        if not line.startswith("AI"):
            line_columns = line.split('\t')
            if float(line_columns[0]) <= 0:
                classification = "N"
            elif float(line_columns[0]) < 30:
                classification = "L"
            else:
                classification = "P"
            gene_classification[line_columns[2]] = classification

    gene_location = {}
    for line in ann_file:
        gene_location[line[3]] = line[0] + "\t" + line[1] + "\t" + line[2] + "\n"

    gene_to_test = []
    tree_results_in = open(tree_results_file, 'r')
    orig_tree_lines = {}

    for line in tree_results_in:
        line = line.rstrip('\n')
        line_columns = line.split('\t')
        gene_to_test.append(line_columns[2])
        orig_tree_lines[line_columns[2]] = line
        if line_columns[0] == "HGT-NT" or line_columns[0] == "HGT":
            gene_classification[line_columns[2]] = "H"
        elif line_columns[0] == "NO" or line_columns[0] == "NO-OT":
            gene_classification[line_columns[2]] = "N"
        elif line_columns[0] == "COMPLEX":
            gene_classification[line_columns[2]] = "C"
        else:
            gene_classification[line_columns[2]] = "-"

    t_prox_file_path = tree_results_file + ".prox.txt"
    t_prox_file = open(t_prox_file_path, "w")

    for gene in gene_to_test:
        gene_loc = BedTool(gene_location[gene], from_string=True)
        nearby_upstream =  gene_loc.closest(ann_file ,id=True, fu=True, io=True, D="a", k=k_value)
        nearby_downstream = gene_loc.closest(ann_file ,iu=True, fd=True, io=True, D="a", k=k_value)
        upstream_string = get_classification_string(gene_classification, nearby_upstream, True)
        downstream_string = get_classification_string(gene_classification, nearby_downstream, False)
        close_string = upstream_string + downstream_string
        score = contamination_score(close_string, k_value)
        t_prox_file.write(orig_tree_lines[gene] + "\t" + upstream_string + "|X|" +downstream_string + "\t" + str(score) + "\n")

def gff_2_bed(gff_file,mode):
    gff_file_in = open(gff_file,'r')
    bed_file = ""
    for line in gff_file_in:
        line = line.rstrip('\n')
        line_columns = line.split('\t')
        if not line.startswith("#") and line_columns[2] == "mRNA" and mode == "Eugene":
            gene = line_columns[8].replace(":", ";").split(";")[1]
            bed_file += line_columns[0] + "\t" + line_columns[3] + "\t" + line_columns[4] + "\t" + gene + "\n"
        if not line.startswith("#") and line_columns[2] == "transcript" and mode == "Augustus":
            gene = line_columns[8]
            bed_file += line_columns[0] + "\t" + line_columns[3] + "\t" + line_columns[4] + "\t" + gene + "\n"
    return bed_file

def get_classification_string(gene_classification,nearby_genes,upstream):
    result_string=""
    for line in nearby_genes:
        if line[6] in gene_classification:
            result_string += gene_classification[line[6]]
        elif line[6] == ".":
            result_string += ""
        else:
            result_string += "-"
    if upstream:
        result_string = result_string[::-1]
    return result_string

def contamination_score(close_string,k_value):
    N_score = (1/(2*k_value))
    H_score = -(1/(2*k_value))
    C_score = 0.0
    U_score = 0.1/k_value
    P_score = -(0.2/k_value)
    L_score = -(0.1/k_value)

    N_count = close_string.count('N')*N_score
    C_count = close_string.count('C')*C_score
    L_count = close_string.count('L')*L_score
    H_count = close_string.count('H')*H_score
    P_count = close_string.count('P')*P_score
    U_count = close_string.count('-')*U_score

    score = N_count + C_count + L_count + H_count + P_count + U_count
    score = round(score,2)
    return(score)

if __name__ == '__main__':
    main()
