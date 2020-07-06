#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""create_recall_plots.py
    Usage: create_recall_plots.py -t <FILE> [-a <INT>]

    Options:
        -h, --help                  show this
        -t, --tree_results <FILE>   tree results file
        -a, --ai_max <INT>          maximum ai for last column [default: 31]
"""


from docopt import docopt
import math
import matplotlib.pyplot as plt


def main():

    args = docopt(__doc__)
    tree_results_file = args['--tree_results']
    ai_max = int(args['--ai_max'])

    y_data = list(range(ai_max))
    list_hgts_wtoi = [0] * ai_max
    list_hgts_all = [0] * ai_max
    list_others = [0] * ai_max
    list_per_complex = [0] * ai_max
    list_per_no = [0] * ai_max

    types = {}
    ais = {}
    num_hgts_wtoi = 0
    num_hgts_all = 0
    num_others = 0
    t_file = open(tree_results_file,'r')
    for line in t_file:
        line = line.rstrip('\n')
        type, tree, gene, ai = line.split('\t')
        ai = float(ai)
        types[gene]=type
        if ai > ai_max:
            ai = ai_max
        ais[gene]=ai

        if type == "HGT":
            num_hgts_wtoi += 1
            num_hgts_all += 1
            for i in range(0,math.ceil(ai),1):
                list_hgts_wtoi[i]+=1
                list_hgts_all[i]+=1

        elif type == "HGT-NT":
            num_hgts_all += 1
            for i in range(0,math.ceil(ai),1):
                list_hgts_all[i]+=1

        elif type == "COMPLEX":
            for i in range(0,math.ceil(ai),1):
                list_others[i]+=1
                list_per_complex[i]+=1
        else:
            for i in range(0,math.ceil(ai),1):
                list_others[i]+=1
                list_per_no[i]+=1
    # For wtoi
    sen_hgts_wtoi = [(x * 100) / num_hgts_wtoi for x in list_hgts_wtoi]
    spe_hgts_wtoi_cons = [ x + z for x,z in zip(list_hgts_wtoi,list_others)]
    spe_hgts_wtoi = [(x * 100) / z for x ,z in zip(list_hgts_wtoi,spe_hgts_wtoi_cons)]

    # For all
    sen_hgts_all = [(x * 100) / num_hgts_all for x in list_hgts_all]
    spe_hgts_all_cons = [ x + z for x,z in zip(list_hgts_all,list_others)]
    spe_hgts_all = [(x * 100) / z for x ,z in zip(list_hgts_all,spe_hgts_all_cons)]

    #Per complex
    per_complex = [(x * 100) / z for x ,z in zip(list_per_complex,spe_hgts_all_cons)]
    per_no = [(x * 100) / z for x ,z in zip(list_per_no,spe_hgts_all_cons)]

    out_fig = tree_results_file + ".pdf"
    plt.plot(sen_hgts_wtoi,label="Sensitivity w TOI")
    plt.plot(spe_hgts_wtoi,label="Specificity w TOI")
    plt.plot(sen_hgts_all,label="Sensitivity All")
    plt.plot(spe_hgts_all,label="Specificity All")
    plt.plot(per_complex,label="Percentage Complex")
    plt.plot(per_no,label="Percentage No HGTs")
    plt.xlabel("AI index")
    plt.ylabel("Percentage")
    plt.legend()
    plt.savefig(out_fig)
    plt.clf()
    plt.cla()
    plt.close()



if __name__ == '__main__':
    main()
