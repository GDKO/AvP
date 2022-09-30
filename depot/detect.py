#!/usr/bin/env python3

"""
  Usage:
    avp detect -i <DIR> -o <DIR> -g <FILE> -c <FILE> [-t <FILE>]
    avp detect -i <DIR> -o <DIR> -g <FILE> -c <FILE> --trees <STR> --trees_suffix <STR> [-t <FILE>]

  Options:
    -h, --help                  show this
    -i, --input <DIR>           folder with aligned fasta files
    -o, --output <DIR>          creates a directory for all output files
    -g, --groups <FILE>         file containing the groups
    -c, --config_file <FILE>    config file
    -t, --taxonomy <FILE>       gid-taxon file for adding lineage to nexus
    --trees <STR>               use pre-calulated trees stored in output_dir/STR
    --trees_suffix <STR>        suffix of the tree files
"""

import sys
import os
import subprocess
from docopt import docopt
from multiprocessing import Pool
from operator import itemgetter
import yaml

from depot.PetIO import check_programs, check_indir, get_outdir, progress, fix_iqtree, run_fasttree, run_iqtree, get_num_jobs
from depot.PetNexus import make_nexus_file, analyze_tree

# Subroutines
def main():

    # =============== #
    # 	PARAMETERS    #
    # =============== #

    args = docopt(__doc__)

    input_dir = check_indir(args['--input'])
    output_dir = get_outdir(args['--output'])
    groups_file = args['--groups']
    config_yaml = args['--config_file']
    taxonomy_file = args['--taxonomy']
    trees = args['--trees']
    trees_suffix = args['--trees_suffix']

    # =============== #
    # 	MAIN          #
    # =============== #

    """
    0. Setting up
    """

    print ("[+] Setting up")

    stream = open(config_yaml, 'r')
    config_opts = yaml.safe_load(stream)
    stream.close()

    data_type = config_opts["data_type"]
    threads = config_opts["max_threads"]
    fastml = config_opts["fastml"]
    node_support = config_opts["node_support"]
    complex_per = [config_opts["complex_per_ingroup"],config_opts["complex_per_donor"]]

    #Check programs
    if fastml:
        check_programs("fasttree")
    else:
        check_programs("iqtree")

    # IQ-TREE parameters
    iq_threads = config_opts["iq_threads"]
    iqmodel = config_opts["iqmodel"]
    ufbootstrap = config_opts["ufbootstrap"]

    fasttree_model = ""

    if data_type == "AA":
        fasttree_model = "-gamma -lg"
    elif data_type == "DNA":
        fasttree_model = "-gamma -gtr"
        iqmodel = "-m GTR+G"
    else:
        print("Throw error\n")

    jobs = get_num_jobs(fastml, threads, iq_threads)

    # Get colors
    colors = {}
    colors_path = os.path.join(sys.path[0],"depot","colors.txt")
    with open(colors_path,'r',encoding='utf8') as fhr_colors:
        for line in fhr_colors:
            if "#" not in line:
                (color,name) = line.rstrip('\n').split(':')
                colors[name] = color

    # Get the lineage
    lineage = {}
    if taxonomy_file:
        f_taxonomy = open(taxonomy_file,'r')

        for line in f_taxonomy:
            line = line.rstrip('\n')
            line_columns = line.split('\t')
            lineage[line_columns[0]] = line_columns[1].split(';')

        f_taxonomy.close()

    # Get the groups
    groups = {}

    f_groups = open(groups_file,'r')

    group_number = 0
    gene_number = 0

    group_dict = {}
    for line in f_groups:
        line = line.rstrip('\n')
        line_columns = line.split('\t')
        if not os.path.isfile(os.path.join(input_dir,line_columns[0])):
            print("\t[!] FATAL ERROR: '{}' not found in fasta directory".format(line_columns[0]))
            sys.exit()
        group_dict[line_columns[0]] = int(line_columns[1])
        groups[line_columns[0]] = line_columns[2:]
        group_number += 1
        gene_number += len(line_columns[2:])

    print ("[!] Found " + str(group_number) + " groups and " + str(gene_number) + " genes")
    f_groups.close()

    """
    1. Trees
    -> FastTree tree (file XXXXX.fasttree)
    -> IqTree tree   (file XXXXX.treefile)
    """
    p = Pool(jobs)

    if trees:
        print ("[!] Skipping Tree Building...")
        tree_path = get_outdir(trees)
    else:
        if fastml:
            print ("[+] Reconstructing phylogenies with FastTree")
            tree_path = get_outdir(output_dir, add_dir="fasttree")
        else:
            print ("[+] Reconstructing phylogenies with IQ-TREE")
            tree_path = get_outdir(output_dir, add_dir="iqtree")

        FNULL = open(os.devnull, 'w')

        i = 0
        job_list = []
        fix_iqtree_list = []

        for group_name, value in sorted(group_dict.items(), key = itemgetter(1), reverse = True):
            group_file = os.path.join(input_dir,group_name)
            if fastml:
                fasttree_params = fasttree_model + " " + group_file + " > " + os.path.join(tree_path, group_name + ".fasttree")
                t_list=[fasttree_params]
                job_list.append(t_list)
            else:
                fname = os.path.join(tree_path,group_name)
                iqtree_params = "-quiet -nt " + str(iq_threads) + " " + str(iqmodel) + " -bb " + str(ufbootstrap) + " -s " + group_file + " -pre " + fname
                job_list.append(iqtree_params)
                fix_iqtree_list.append(fname + ".treefile")

        if fastml:
            for i, _ in enumerate(p.imap_unordered(run_fasttree,job_list),1):
                progress(i,1,len(job_list))
        else:
            for i, _ in enumerate(p.imap_unordered(run_iqtree,job_list),1):
                progress(i,1,len(job_list))
            #Add @ back to the treefile because iqtree removes it
            for tree in fix_iqtree_list:
                fix_iqtree(tree)

    """
    2. Analyzing tree results
    """
    if trees:
        print ("[+] Analyzing " + trees + " results")
        tree_nexus_path = get_outdir(output_dir, add_dir=trees_suffix+"_nexus")
        t_res_path = os.path.join(output_dir, trees_suffix+"_tree_results.txt")
        g_res_path = os.path.join(output_dir, trees_suffix+"_general_results.txt")
    elif fastml:
        print ("[+] Analyzing fasttree results")
        tree_nexus_path = get_outdir(output_dir, add_dir="fasttree_nexus")
        t_res_path = os.path.join(output_dir, "fasttree_tree_results.txt")
        g_res_path = os.path.join(output_dir, "fasttree_general_results.txt")
    else:
        print ("[+] Analyzing IQ-TREE results")
        tree_nexus_path = get_outdir(output_dir, add_dir="iqtree_nexus")
        t_res_path = os.path.join(output_dir, "iqtree_tree_results.txt")
        g_res_path = os.path.join(output_dir, "iqtree_general_results.txt")

    t_res = open(t_res_path,'w')

    result_nohgt = "no_hgt"
    result_unknown = "unknown_topology"
    result_complex = "complex_topology"
    result_hgt = "hgt"

    no_hgt_count = 0
    unknown_count = 0
    complex_count = 0
    strong_hgt_count = 0
    no_ingroup_count = 0
    only_ingroup_count = 0

    i = 0
    for group in groups.keys():
        group_file = os.path.join(input_dir,group)

        if trees:
            phylogeny_file = os.path.join(tree_path,group + "." + trees_suffix)
        elif fastml:
            phylogeny_file = os.path.join(tree_path,group + ".fasttree")
        else:
            phylogeny_file = os.path.join(tree_path,group + ".treefile")

        for gene in groups[group]:
            gene_nexus_path = os.path.join(tree_nexus_path,gene + ".nexus")
            make_nexus_file(gene, group, lineage, gene_nexus_path, group_file, phylogeny_file, colors)
            result_tree = analyze_tree(phylogeny_file,gene+"@StudiedOrganism", node_support, complex_per)
            if (result_tree == result_nohgt):
                no_hgt_count += 1
                t_res.write("NO\t" + group_file + "\t" + phylogeny_file + "\t" + gene + "\n")
            elif (result_tree == "only_Ingroup"):
                no_hgt_count += 1
                only_ingroup_count += 1
                t_res.write("NO-OT\t" + group_file + "\t" + phylogeny_file + "\t" + gene + "\n")
            elif (result_tree == result_unknown):
                unknown_count += 1
                t_res.write("UNKNOWN\t" + group_file + "\t" + phylogeny_file + "\t" + gene + "\n")
            elif (result_tree == result_complex):
                complex_count += 1
                t_res.write("COMPLEX\t" + group_file + "\t" + phylogeny_file + "\t" + gene + "\n")
            elif (result_tree == result_hgt):
                strong_hgt_count += 1
                t_res.write("HGT\t" + group_file + "\t" + phylogeny_file + "\t" + gene + "\n")
            elif (result_tree == "no_Ingroup"):
                strong_hgt_count += 1
                no_ingroup_count += 1
                t_res.write("HGT-NT\t" + group_file + "\t" + phylogeny_file + "\t" + gene + "\n")
            else:
                print("\t[!] FATAL ERROR: Code 1")
                sys.exit()

            i = i + 1
            progress(i,1,gene_number)

    t_res.close()

    """
    X. Summary file
    """

    # Final output file
    g_res = open(g_res_path,'w')
    g_res.write("Proteins analyzed\t: " + str(gene_number) + "\n")
    g_res.write("Proteins with no Ingroup\t: " + str(no_ingroup_count) + "\n")
    g_res.write("Proteins with only Ingroup\t: " + str(only_ingroup_count) + "\n")
    g_res.write("\n")
    g_res.write("Unknown Topology\t: " + str(unknown_count) + "\n")
    g_res.write("No HGT support\t\t: " + str(no_hgt_count) + "\n")
    g_res.write("Complex topology\t: " + str(complex_count) + "\n")
    g_res.write("Strong HGT support\t: " + str(strong_hgt_count) + "\n")
    g_res.close()

if __name__ == '__main__':
    main()
