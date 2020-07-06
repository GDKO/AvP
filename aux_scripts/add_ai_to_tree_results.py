#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""add_ai_to_tree_results.py
    Usage: add_ai_to_tree_results.py -a <FILE> -t <FILE>

    Options:
        -h, --help                  show this
        -a, --aifeatures <FILE>     alienness features file
        -t, --tree_results <FILE>   tree results file
"""

import os
from docopt import docopt


def main():

    args = docopt(__doc__)
    ai_features = args['--aifeatures']
    tree_results_file = args['--tree_results']

    ai = {}

    ai_file = open(ai_features,'r')
    for line in ai_file:
        line = line = line.rstrip('\n')
        line_columns = line.split('\t')
        ai[line_columns[2]]=line_columns[0]

    t_file = open(tree_results_file,'r')
    t_ai_file_path = tree_results_file + ".ai.txt"
    t_ai_file = open(t_ai_file_path,'w')
    for line in t_file:
        line = line.rstrip('\n')
        line_columns = line.split('\t')
        t_ai_file.write(line + '\t' + str(ai[line_columns[2]]) + '\n')

    t_ai_file.close()

if __name__ == '__main__':
    main()
