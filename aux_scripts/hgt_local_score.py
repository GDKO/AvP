#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""closest_genes.py
    Usage: closest_genes.py -f <FILE> -a <FILE> -t <FILE> -m <INT> [-k <INT>]

    Options:
        -h, --help                  show this
        -f, --file <FILE>           gff or bed file with gene coordinates
        -a, --ai_file <FILE>        Alienness index full file
        -t, --tree_results <FILE>   Tree results file
        -k, --k_value <INT>         number of close genes to process [default: 5]
        -m, --mode <INT>            mode gff, check documentation

"""

import os
import sys
import re
from docopt import docopt

depot_path = os.path.join(sys.path[0], "../")
sys.path.insert(0,depot_path)
from depot.PetIO import open_file

def main():
    args = docopt(__doc__)
    genome_annotation_file = args['--file']
    ai_file = args['--ai_file']
    tree_results_file = args['--tree_results']
    k_value = int(args['--k_value'])
    mode = int(args['--mode'])

    # mode 0
    # [scaffold]  [source]  mRNA  [start]  [end]  [score]  [strand]  [frame]  ID=[transcript];Parent=[gene]
    mode_gff = [["mRNA",";","="]]

    scaf_location, gene_location = gff3_to_dict(genome_annotation_file,mode_gff[mode])
    ai_file_in = open(ai_file, 'r')
    gene_classification = {}
    for line in ai_file_in:
        if not line.startswith("query"):
            line_columns = line.split('\t')
            gene = line_columns[0]
            ai = line_columns[3]
            if float(ai) <= 0:
                classification = "N"
            elif float(ai) < 30:
                classification = "L"
            else:
                classification = "P"
            gene_classification[gene] = classification


    gene_to_test = []
    tree_results_in = open(tree_results_file, 'r')
    orig_tree_lines = {}

    for line in tree_results_in:
        line = line.rstrip('\n')
        line_columns = line.split('\t')
        gene = line_columns[3]
        orig_tree_lines[gene] = line
        if line_columns[0] == "HGT-NT" or line_columns[0] == "HGT":
            gene_classification[gene] = "H"
            gene_to_test.append(gene)
        elif line_columns[0] == "NO" or line_columns[0] == "NO-OT":
            gene_classification[gene] = "N"
        elif line_columns[0] == "COMPLEX":
            gene_classification[gene] = "C"
        else:
            gene_classification[gene] = "-"

    t_prox_file_path = tree_results_file + ".prox.txt"
    t_prox_file = open(t_prox_file_path, "w")

    for gene in gene_to_test:
        gene_scaf = gene_location[gene]
        index = int(scaf_location[gene_scaf].index(gene))
        nearby_upstream = scaf_location[gene_scaf][max(0,index-k_value):index]
        nearby_downstream = scaf_location[gene_scaf][index+1:min(len(scaf_location[gene_scaf]),index+k_value)+1]
        upstream_string = get_classification_string(gene_classification, nearby_upstream)
        downstream_string = get_classification_string(gene_classification, nearby_downstream)
        close_string = upstream_string + downstream_string
        score = contamination_score(close_string, k_value)
        t_prox_file.write(orig_tree_lines[gene] + "\t" + upstream_string + "|X|" +downstream_string + "\t" + str(score) + "\n")

def gff3_to_dict(gff_file,mode_gff):
    scaf_location = {}
    sorted_scaf_location = {}
    gene_location = {}
    with open_file(gff_file) as fhr_gff:
        for line in fhr_gff:
            if not line.startswith("#"):
                gff_column = line.rstrip('\n').split('\t')
                scaf = gff_column[0]
                feature = gff_column[2]
                start = gff_column[3]
                attr = gff_column[8]
                transcript = attr.split(mode_gff[1])[0].split(mode_gff[2])[1]
                if feature == mode_gff[0]:
                    if not scaf in scaf_location.keys():
                        scaf_location[scaf] = {}
                        scaf_location[scaf]["start"] = []
                        scaf_location[scaf]["attr"] = []
                        sorted_scaf_location[scaf] = {}
                    scaf_location[scaf]["start"].append(int(start))
                    scaf_location[scaf]["attr"].append(transcript)
                    gene_location[transcript] = scaf

    for scaf in scaf_location:
        X = scaf_location[scaf]["start"]
        Y = scaf_location[scaf]["attr"]
        L = [x for _,x in sorted(zip(X,Y))]
        sorted_scaf_location[scaf] = L
    return sorted_scaf_location, gene_location


def get_classification_string(gene_classification,nearby_genes):
    result_string = ""
    for gene in nearby_genes:
        if gene in gene_classification:
            result_string += gene_classification[gene]
        else:
            result_string += "-"
    return result_string

def contamination_score(close_string,k_value):
    prox_count = len(close_string)
    f_value = 1
    if prox_count != 0:
        n=(prox_count-(2*k_value))/(2*k_value)
        f_value = math.exp(n)
    N_score = (1/(2*k_value*f_value))
    H_score = -(1/(2*k_value*f_value))
    C_score = 0.0
    U_score = 0.1/(k_value*f_value)
    P_score = -(0.2/(k_value*f_value))
    L_score = -(0.1/(k_value*f_value))

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
