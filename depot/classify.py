#!/usr/bin/env python3

"""
  Usage:
    avp classify -i <DIR> -t <FILE> -f <FILE> -c <FILE> -o <DIR>

  Options:
    -h, --help                    show this
    -i, --input_nexus_dir <DIR>   folder with nexus files
    -t, --tree_results <FILE>     tree results file
    -f, --classification <FILE>   classification file
    -c, --config_file <FILE>      config file
    -o, --output <DIR>            creates a directory for all output files
"""

import shutil
import os
import yaml
from docopt import docopt
from depot.PetIO import get_outdir, progress, check_indir
from depot.PetNexus import classify_tree


def main():

    # =============== #
    # 	PARAMETERS    #
    # =============== #

    args = docopt(__doc__)

    input_nexus_dir = check_indir(args['--input_nexus_dir'])
    tree_results_file = args['--tree_results']
    classification_file = args['--classification']
    config_yaml = args['--config_file']
    output_dir = get_outdir(args['--output'])

    out_path = get_outdir(output_dir, add_dir="classification")

    stream = open(config_yaml, 'r')
    config_opts = yaml.safe_load(stream)
    stream.close()

    node_support = config_opts["node_support"]
    complex_per_node = config_opts["complex_per_node"]

    print ("[+] Setting up")

    classification = {}
    results = {}
    c_file = open(classification_file,'r')

    for line in c_file:
        if not line.startswith("#"):
            line = line.rstrip('\n')
            line_columns = line.split('\t')
            classification[line_columns[0]] = line_columns[1].split(';')

    if len(classification) > 1 :
        hgt_complex_dir = get_outdir(out_path, add_dir="HGT_complex")
        results["HGT_complex"] = 0

    full_ranks = []
    for rank in classification:
        if len(classification[rank]) > 1:
            get_outdir(out_path, add_dir= rank + "_complex")
            results[rank + "_complex"] = 0
        for inside_rank in classification[rank]:
            get_outdir(out_path, add_dir=inside_rank)
            results[inside_rank] = 0
            full_ranks.append(inside_rank)

    results["Complex"] = 0
    results["Ingroup"] = 0
    results["Unknown"] = 0
    complex_dir = get_outdir(out_path, add_dir= "Complex")
    ingroup_dir = get_outdir(out_path, add_dir= "Ingroup")
    unknown_dir = get_outdir(out_path, add_dir= "Unknown")

    tree = {}
    type = {}
    rank_results = {}

    t_file = open(tree_results_file,'r')

    for line in t_file:
        line = line.rstrip('\n')
        line_columns = line.split('\t')
        type[line_columns[3]] = line_columns[0]
        tree[line_columns[3]] = line_columns[2]

    print ("[+] Classifying Trees")
    i = 0
    for gene in tree:
        if type[gene] == "COMPLEX":
            shutil.copy(input_nexus_dir+gene+".nexus",complex_dir)
            results["Complex"] +=1
            rank_results[gene] = "Complex"
        elif type[gene] == "NO" or type[gene] == "NO-OT":
            shutil.copy(input_nexus_dir+gene+".nexus",ingroup_dir)
            results["Ingroup"] += 1
            rank_results[gene] = "Ingroup"
        elif type[gene] == "UNKNOWN":
            shutil.copy(input_nexus_dir+gene+".nexus",unknown_dir)
            results["Unknown"] += 1
            rank_results[gene] = "Unknown"
        else:
            final_nodes = classify_tree(tree[gene],gene+"@StudiedOrganism",full_ranks, complex_per_node)
            if len(final_nodes) == 0:
                shutil.copy(input_nexus_dir+gene+".nexus",unknown_dir)
                results["Unknown"] += 1
                rank_results[gene] = "Unknown"
            elif len(final_nodes) == 1:
                shutil.copy(input_nexus_dir+gene+".nexus",os.path.join(out_path,final_nodes[0]))
                results[final_nodes[0]] += 1
                rank_results[gene] = final_nodes[0]
            else:
                final_ranks = []
                for node in final_nodes:
                    for rank in classification:
                        if node in classification[rank]:
                            final_ranks.append(rank)
                if len(list(set(final_ranks))) == 0:
                    shutil.copy(input_nexus_dir+gene+".nexus",unknown_dir)
                    results["Unknown"] += 1
                    rank_results[gene] = "Unknown"
                elif len(list(set(final_ranks))) == 1:
                    shutil.copy(input_nexus_dir+gene+".nexus",os.path.join(out_path,final_ranks[0]+"_complex"))
                    results[final_ranks[0]+"_complex"] += 1
                    rank_results[gene] = final_ranks[0]+"_complex"
                else:
                    shutil.copy(input_nexus_dir+gene+".nexus",hgt_complex_dir)
                    results["HGT_complex"] += 1
                    rank_results[gene] = "HGT_complex"
        i = i + 1
        progress(i,1,len(tree))

    """
    X. Summary file
    """

    t_res_path = os.path.join(output_dir, "classification_tree_results.txt")
    t_res = open(t_res_path,'w')
    for gene in tree:
        t_res.write(rank_results[gene]+'\t'+tree[gene]+'\t'+gene+'\n')

    g_res_path = os.path.join(output_dir, "classification_results.txt")
    g_res = open(g_res_path,'w')
    for item in results:
        l = int(len(item)/8)
        if l == 0:
            tabs = '\t\t\t'
        elif l == 1:
            tabs = '\t\t'
        else:
            tabs = '\t'
        g_res.write(item+tabs+": "+str(results[item])+'\n')

    """
    Clean empty directories
    """
    #Clean empty direcotries
    for dir in os.listdir(out_path):
        try:
            os.rmdir(os.path.join(out_path,dir))
        except OSError as ex:
            pass


if __name__ == '__main__':
    main()
