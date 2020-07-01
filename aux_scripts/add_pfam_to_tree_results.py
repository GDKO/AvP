#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""add_pfam_to_tree_results.py
    Usage: add_pfam_to_tree_results.py -p <FILE> -t <FILE> [-e <STR>]

    Options:
        -h, --help                  show this
        -p, --pfam_results <FILE>   pfam results file from hmmer3
        -t, --tree_results <FILE>   tree results file
        -e, --evalue_cutoff <STR>   evalue_cutoff [default: 1e-5]
        --version                   print version
"""

import os
from docopt import docopt


def main():

    args = docopt(__doc__,version='1.00')
    pfam_file = args['--pfam_results']
    tree_results_file = args['--tree_results']
    evalue_cutoff = args['--evalue_cutoff']

    rename_prot = True
    str_wrong = ":"
    str_good = "_"    

    pfam = {}

    pf_file = open(pfam_file,'r')
    for line in pf_file:
        line = line = line.rstrip('\n')
        if not line.startswith("#"):
            line_columns = line.split(None,maxsplit=22)
            if rename_prot:
                transcript_name = line_columns[3].replace(str_wrong,str_good)
            else:
                transcript_name = line_columns[3]
            if transcript_name not in pfam:
                pfam[transcript_name]=[]
            if line_columns[6]<evalue_cutoff:
                pfam[transcript_name].append(line_columns[22])


    t_file = open(tree_results_file,'r')
    t_pf_file_path = tree_results_file + ".pfam.txt"
    t_pf_file = open(t_pf_file_path,'w')
    for line in t_file:
        line = line = line.rstrip('\n')
        line_columns = line.split('\t')
        transcript_name = line_columns[2]
        if transcript_name in pfam:
            pfam_list = list(set(pfam[transcript_name]))
            pfam_list_string = '\t'.join(pfam_list)
        else:
            pfam_list_string = ""
        t_pf_file.write(line + '\t' + str(pfam_list_string) + '\n')

    t_pf_file.close()

if __name__ == '__main__':
    main()
