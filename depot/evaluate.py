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
import yaml
import subprocess

from docopt import docopt
from Bio import SeqIO
from multiprocessing import Pool
from depot.PetIO import check_programs, open_file, get_outdir, progress, check_indir, run_fasttree, run_iqtree, fix_iqtree, get_num_jobs

def main():

    # =============== #
    # 	PARAMETERS    #
    # =============== #

    args = docopt(__doc__)

    input_dir = check_indir(args['--input'])
    tree_results_file = args['--tree_results']
    output_dir = get_outdir(args['--output'])
    config_yaml = args['--config_file']

    stream = open(config_yaml, 'r')
    config_opts = yaml.safe_load(stream)
    stream.close()

    data_type = config_opts["data_type"]
    threads = config_opts["max_threads"]
    fastml = config_opts["fastml"]

    iq_threads = config_opts["iq_threads"]
    iqmodel = config_opts["iqmodel"]
    ufbootstrap = config_opts["ufbootstrap"]

    fasttree_model = ""
    fastml_iqmodel = ""

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

    p = Pool(jobs)

    # Ingroup, Gene, Exclude [1]
    # Unknown, Rest of Genes [-] means freely
    # Rest [0]

    t_file = open(tree_results_file,'r')

    tree = {}
    type = {}
    alignment = {}

    for line in t_file:
        line = line.rstrip('\n')
        line_columns = line.split('\t')
        type[line_columns[3]] = line_columns[0]
        alignment[line_columns[3]] = line_columns[1]
        tree[line_columns[3]] = line_columns[2]

    fasttree_job_list = []
    iqtree_ml_job_list = []
    fix_iqtree_list = []

    con_tree_list = []
    iqtree_alt_job_list = []

    selected_genes = []

    if fastml:
        print ("[+] Reconstructing constrainted topologies with FastTree")
    else:
        print ("[+] Reconstructing constrainted topologies with IQ-TREE")

    for gene in tree:
        if str(type[gene]) == "HGT":

            selected_genes.append(gene)

            ingroup = ["EGP", "Ingroup"]
            ingroup.append(gene)
            filename = alignment[gene]

            ingroup_list = []
            free_list = []
            rest_list = []

            with open_file(filename) as handle:
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
            con = open(con_path,'w')

            FNULL = open(os.devnull, 'w')
            all_trees_path =  os.path.join(out_path,gene + ".trees")

            if fastml:
                # Constaint file for FastTree
                total = len(ingroup_list) + len(free_list) + len(rest_list)
                con.write(" " + str(total) +" 1\n")
                for gene_l in ingroup_list:
                    con.write(gene_l + " 1\n")
                for gene_l in free_list:
                    con.write(gene_l + " -\n")
                for gene_l in rest_list:
                    con.write(gene_l + " 0\n")
                con.close()
                # FastTree params
                tree_path = os.path.join(out_path,gene + ".fasttree")
                fasttree_params = fasttree_model + " -constraints " + con_path + " " + alignment[gene] + " > " + tree_path
                t_list=[fasttree_params]
                fasttree_job_list.append(t_list)

                con_tree_list.append([tree_path,tree[gene],all_trees_path])
                #IqTree Alt Params
                iqtree_alt_params = "-quiet -nt " + str(iq_threads) + " " + fastml_iqmodel + " -n 0 -zb 10000 -au -pre " + tree_path + " -s " + alignment[gene] + " -z " + all_trees_path
                iqtree_alt_job_list.append(iqtree_alt_params)

            else:
                # Constaint file for IqTree
                ingroup_list = [w.replace('@', '_') for w in ingroup_list]
                rest_list = [w.replace('@', '_') for w in rest_list]
                con.write("((" + ', '.join(ingroup_list) + ")," + "" + ', '.join(rest_list) + ");" + "\n")
                con.close()
                #IqTree Params
                tree_pre = os.path.join(out_path,gene)
                tree_path = tree_pre + ".treefile"
                iqtree_ml_params = "-quiet -nt " + str(iq_threads) + " " + iqmodel + " -bb " + str(ufbootstrap) + " -g " + con_path + " -s " + alignment[gene] + " -pre " + tree_pre
                iqtree_ml_job_list.append(iqtree_ml_params)
                fix_iqtree_list.append(tree_path)

                con_tree_list.append([tree_path,tree[gene],all_trees_path])
                #IqTree Alt Params
                iqtree_alt_params = "-quiet -nt " + str(iq_threads) + " " + iqmodel + " -n 0 -zb 10000 -au -pre " + tree_path + " -s " + alignment[gene] + " -z " + all_trees_path
                iqtree_alt_job_list.append(iqtree_alt_params)

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


    #Add jobs for iqtree
    jobs = get_num_jobs(False, threads, iq_threads)
    p = Pool(jobs)

    print ("[+] Calculating stats with IQ-TREE")

    for i, _ in enumerate(p.imap_unordered(run_iqtree,iqtree_alt_job_list),1):
        progress(i,1,len(iqtree_alt_job_list))

    alt_res_path = os.path.join(output_dir, "alt_topology_results.tsv")
    alt_res = open(alt_res_path,'w')
    alt_res.write("#gene"+ "\t" +"deltaL" + "\t" + "p-AU"+ "\t" + "Significantly worse" + "\n")
    #Create final file
    for gene in selected_genes:

        if fastml:
            file = out_path + "/" + gene + ".fasttree.iqtree"
        else:
            file = out_path + "/" + gene + ".treefile.iqtree"

        with open(file) as f:
            content = f.readlines()

        for i, line in enumerate(content):
            if "USER TREES" in line:
                tree_1 = content[i+7].split()
                tree_2 = content[i+8].split()
        deltaL = float(tree_1[2])-float(tree_2[2])
        au=tree_1[11]
        if tree_1[12] == "-":
            sf=1
        else:
            sf =0
        alt_res.write(gene + "\t" + str(deltaL) + "\t" + str(au) + "\t" + str(sf)+"\n")

    alt_res.close()



def concatenate_trees(con_tree_list):
    tree_path, tree_gene, all_trees_path = con_tree_list
    FNULL = open(os.devnull, 'w')
    subprocess.call("cat " + tree_path + " " + tree_gene + " | sed 's/@/_/g' > " + all_trees_path, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

if __name__ == '__main__':
    main()
