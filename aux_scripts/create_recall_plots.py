#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""create_recall_plots.py
    Usage: create_recall_plots.py -t <FILE> [-a <INT>]

    Options:
        -h, --help                  show this
        -t, --tree_results <FILE>   tree results file
        -a, --ai_max <INT>          maximum ai for last column [default: 40]
"""


from docopt import docopt
import math
import matplotlib.pyplot as plt
import numpy as np


def main():

    args = docopt(__doc__)
    tree_results_file = args['--tree_results']
    ai_max = int(args['--ai_max'])

    t_file = open(tree_results_file,'r')
    ai_list = []
    hgt_list = []
    for line in t_file:
        line = line.rstrip('\n')
        line_columns = line.split('\t')
        ai = float(line_columns[4])
        ai_list.append(ai)
        if "HGT" in line_columns[0]: # and "HGT-NT" not in line_columns[0]:
            hgt_list.append(ai)

    ai_group_list = []
    hgt_group_list = []
    for ai in range(0,ai_max+1):
        ai_group_list.append(sum(i >= ai for i in ai_list))
        hgt_group_list.append(sum(i >= ai for i in hgt_list))

    f1_list = []
    ppv_list = []
    sen_list = []
    diff_list = []
    for ai in range(0,ai_max+1):
        f1_score = 2*hgt_group_list[ai]/(hgt_group_list[0]+ai_group_list[ai])
        ppv_score = hgt_group_list[ai]/ai_group_list[ai]
        sen_score = hgt_group_list[ai]/hgt_group_list[0]
        f1_list.append(f1_score)
        ppv_list.append(ppv_score)
        sen_list.append(sen_score)
        diff_list.append(abs(ppv_score-sen_score))

    max_f1 = f1_list.index(max(f1_list))

    out_fig = tree_results_file + ".pdf"

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ppv_list,label="Precision",color="#785EF0")
    ax.plot(sen_list,label="Sensitivity",color="#DC267F")
    ax.plot(f1_list,label="F1 score",color="#FFB000")
    ax.vlines(max_f1,ppv_list[0],1,color="#FFB000",linestyles="dashed",linewidth=1)
    ax.set(xlabel='AI index')
    plt.xticks(fontsize=6)
    plt.yticks(fontsize=6)
    x0, x1 = ax.get_xlim()
    visible_ticks = [t for t in ax.get_xticks() if t>=x0 and t<=x1]
    x_ticks = np.append(visible_ticks, max_f1)
    ax.set_xticks(x_ticks)
    ax.legend()
    fig.savefig(out_fig)


if __name__ == '__main__':
    main()
