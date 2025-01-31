#!/usr/bin/env python3

"""
  Usage:
    avp evaluate -i <DIR> -t <FILE> -o <DIR> -c <FILE>

  Options:
    -h, --help                    show this
    -i, --input <DIR>             folder with aligned fasta files
    -t, --tree_results <FILE>     tree results file
    -o, --output <DIR>            creates a directory for all output files
    -c, --config_file <FILE>      config file
"""

import sys
import os
from multiprocessing import Pool
import subprocess
import yaml

from docopt import docopt
from Bio import SeqIO
from depot.PetIO import check_programs, open_file, get_outdir, progress, check_indir, run_fasttree, run_iqtree, fix_iqtree, get_num_jobs

def main():
    """Evaluates alternative topologies"""

    # =============== #
    # 	PARAMETERS    #
    # =============== #

    args = docopt(__doc__)

    check_indir(args['--input'])
    tree_results_file = args['--tree_results']
    output_dir = get_outdir(args['--output'])
    config_yaml = args['--config_file']

    with open(config_yaml, 'r', encoding = "UTF-8") as stream:
        config_opts = yaml.safe_load(stream)

    data_type = config_opts["data_type"]
    threads = config_opts["max_threads"]
    fastml = config_opts["fastml"]

    iq_threads = config_opts["iq_threads"]
    iqmodel = config_opts["iqmodel"]
    ufbootstrap = config_opts["ufbootstrap"]

    if data_type == "AA":
        fasttree_model = "-gamma -lg"
        fastml_iqmodel = "-m LG+G"
    elif data_type == "DNA":
        fasttree_model = "-gamma -gtr"
        iqmodel = "-m GTR+G"
        fastml_iqmodel = iqmodel
    else:
        sys.exit("data_type should be either AA or DNA")

    jobs = get_num_jobs(fastml, threads, iq_threads)

    out_path = get_outdir(output_dir, add_dir="alt_topology")

    print ("[+] Setting up")

    if fastml:
        check_programs("fasttree", "iqtree")
    else:
        check_programs("iqtree")


    #tree_results_d = [type_res, alignment_file, tree_gene]
    tree_results_d = parse_tree_results_file(tree_results_file)

    # Ingroup, Gene, Exclude [1]
    # Unknown, Rest of Genes [-] means freely
    # Rest [0]

    con_tree_list = []
    iqtree_alt_job_list = []

    selected_genes = []

    if fastml:
        print ("[+] Reconstructing constrainted topologies with FastTree")
        fasttree_job_list = []
    else:
        print ("[+] Reconstructing constrainted topologies with IQ-TREE")
        iqtree_ml_job_list = []
        fix_iqtree_list = []

    for gene, tree_results_l in tree_results_d.items():
        type_res, alignment_file, tree_gene = tree_results_l
        if type_res == "HGT":

            selected_genes.append(gene)

            ingroup = ["EGP", "Ingroup"]
            ingroup.append(gene)

            ingroup_list = []
            free_list = []
            rest_list = []

            with open_file(alignment_file) as handle:
                for record in SeqIO.parse(handle,"fasta"):
                    if any(n in record.id for n in ingroup):
                        ingroup_list.append(record.id)
                    elif "Unknown" in record.id or "StudiedOrganism" in record.id:
                        free_list.append(record.id)
                    else:
                        rest_list.append(record.id)
            #Set up files
            all_trees_path =  os.path.join(out_path,gene + ".trees")
            con_path = os.path.join(out_path, gene + ".constraint")

            with open(con_path, 'w', encoding = "UTF-8") as con:
                if fastml:
                    # Constaint file for FastTree
                    total = len(ingroup_list) + len(free_list) + len(rest_list)
                    con.write(f" {total} 1\n")
                    for gene_l in ingroup_list:
                        con.write(f"{gene_l} 1\n")
                    for gene_l in free_list:
                        con.write(f"{gene_l} -\n")
                    for gene_l in rest_list:
                        con.write(f"{gene_l} 0\n")
                    # FastTree params
                    tree_path = os.path.join(out_path,gene + ".fasttree")
                    fasttree_params = f"{fasttree_model} -constraints {con_path} {alignment_file} > {tree_path}"
                    t_list=[fasttree_params]
                    fasttree_job_list.append(t_list)

                    con_tree_list.append([tree_path,tree_gene,all_trees_path])
                    #IqTree Alt Params
                    iqtree_alt_params = f"-quiet -nt {iq_threads} {fastml_iqmodel} -n 0 -zb 10000 -au -pre {tree_path} -s {alignment_file} -z {all_trees_path}"
                    iqtree_alt_job_list.append(iqtree_alt_params)

                else:
                    # Constaint file for IqTree
                    ingroup_list = [w.replace('@', '_') for w in ingroup_list]
                    rest_list = [w.replace('@', '_') for w in rest_list]
                    con.write(f"(({', '.join(ingroup_list)}),{', '.join(rest_list)});\n")
                    #IqTree Params
                    tree_pre = os.path.join(out_path,gene)
                    tree_path = tree_pre + ".treefile"
                    iqtree_ml_params = f"-quiet -nt {iq_threads} {iqmodel} -bb {ufbootstrap} -g {con_path} -s {alignment_file} -pre {tree_pre}"
                    iqtree_ml_job_list.append(iqtree_ml_params)
                    fix_iqtree_list.append(tree_path)

                    con_tree_list.append([tree_path,tree_gene,all_trees_path])
                    #IqTree Alt Params
                    iqtree_alt_params = f"-quiet -nt {iq_threads} {iqmodel} -n 0 -zb 10000 -au -pre {tree_path} -s {alignment_file} -z {all_trees_path}"
                    iqtree_alt_job_list.append(iqtree_alt_params)

    with Pool(jobs) as p:
        if fastml:
            for i, _ in enumerate(p.imap_unordered(run_fasttree,fasttree_job_list),1):
                progress(i,1,len(fasttree_job_list))
        else:
            for i, _ in enumerate(p.imap_unordered(run_iqtree,iqtree_ml_job_list),1):
                progress(i,1,len(iqtree_ml_job_list))
            #Add @ back to the treefile because iqtree removes it
            for tree in fix_iqtree_list:
                fix_iqtree(tree)

    for con_tree in con_tree_list:
        concatenate_trees(con_tree)

    print ("[+] Calculating stats with IQ-TREE")
    #Add jobs for iqtree
    jobs = get_num_jobs(False, threads, iq_threads)
    with Pool(jobs) as p:

        for i, _ in enumerate(p.imap_unordered(run_iqtree,iqtree_alt_job_list),1):
            progress(i,1,len(iqtree_alt_job_list))

    alt_res_path = os.path.join(output_dir, "alt_topology_results.tsv")
    write_alt_results(alt_res_path, selected_genes, fastml, out_path)

def concatenate_trees(con_tree_list):
    """Concatenate the two trees together"""
    tree_path, tree_gene, all_trees_path = con_tree_list
    with open(os.devnull, 'w', encoding = "UTF-8") as fn:
        subprocess.call("cat " + tree_path + " " + tree_gene + " | sed 's/@/_/g' > " + all_trees_path, shell=True, stdout=fn, stderr=subprocess.STDOUT)

def parse_tree_results_file(tree_results_file):
    """ Parses the tree_results file """
    tree_results_d = {}

    with open(tree_results_file, 'r', encoding = "UTF-8") as t_file:
        for line in t_file:
            line = line.rstrip('\n')
            line_columns = line.split('\t')
            if not os.path.exists(line_columns[2]):
                sys.exit(f"[x] Cannot find tree {line_columns[2]} from {tree_results_file}")
            gene = line_columns[3]
            l = [line_columns[0], line_columns[1], line_columns[2]]
            tree_results_d[gene] = l

    return tree_results_d

def write_alt_results(alt_res_path, selected_genes, fastml, out_path):
    """ Write the alt results tsv """
    with open(alt_res_path, 'w', encoding = "UTF-8") as alt_res:
        alt_res.write("#gene\tdeltaL\tp-AU\tSignificantly worse\n")
        #Create final file
        for gene in selected_genes:

            if fastml:
                file = out_path + "/" + gene + ".fasttree.iqtree"
            else:
                file = out_path + "/" + gene + ".treefile.iqtree"

            with open(file, 'r', encoding = "UTF-8") as f:
                content = f.readlines()

            tree_1 = ""
            tree_2 = ""
            for i, line in enumerate(content):
                if "USER TREES" in line:
                    tree_1 = content[i+7].split()
                    tree_2 = content[i+8].split()
            delta_l = float(tree_1[2])-float(tree_2[2])
            if tree_1[12] == "-":
                sf = 1
            else:
                sf = 0
            alt_res.write(f"{gene}\t{delta_l}\t{tree_1[11]}\t{sf}\n")

if __name__ == '__main__':
    main()
