#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""add_expression_to_tree_results.py
    Usage: add_expression_to_tree_results.py -e <FILE> -t <FILE>

    Options:
        -h, --help                   show this
        -e, --expression_file <FILE> expression_file file
        -t, --tree_results <FILE>    tree results file
        --version                    print version
"""

import os
from docopt import docopt


def main():

    args = docopt(__doc__,version='0.9.0')
    expression_file = args['--expression_file']
    tree_results_file = args['--tree_results']

    trandecoder = True
    trandecoder_split = "__"

    expr = {}

    expr_file = open(expression_file,'r')
    for line in expr_file:
        if not line.startswith("#"):
            line = line = line.rstrip('\n')
            line_columns = line.split('\t')
            expr[line_columns[0]] = []
            for rep_expr in line_columns[1:]:
                expr[line_columns[0]].append(rep_expr)

    t_file = open(tree_results_file,'r')
    t_expr_file_path = tree_results_file + ".expr.txt"
    t_expr_file = open(t_expr_file_path,'w')
    for line in t_file:
        line = line = line.rstrip('\n')
        line_columns = line.split('\t')
        if trandecoder:
            transcript_name = line_columns[2].split(trandecoder_split)[1]
        else:
            transcript_name = line_columns[2]
        expr_list = '\t'.join(expr[transcript_name])
        t_expr_file.write(line + '\t' + str(expr_list) + '\n')

    t_expr_file.close()

if __name__ == '__main__':
    main()
