#!/usr/bin/env python3

"""test_alt_topology.py
  Usage:
    test_alt_topology.py -i <DIR> -t <FILE> -m <STR> -o <DIR> [-x <INT>]

  Options:
    -h, --help                    show this
    -i, --input <DIR>             folder with aligned fasta files
    -t, --tree_results <FILE>     tree results file
    -m, --mode <STR>              select from fasttree, iqtree
    -o, --output <DIR>            creates a directory for all output files
    -x, --threads <INT>           number of threads [default: 2]
    --version                     print version
"""

import sys
import os
depot_path = os.path.join(sys.path[0], "../")
sys.path.insert(0,depot_path)

import subprocess

from docopt import docopt
from Bio import SeqIO
from multiprocessing import Pool
from depot.PetIO import check_programs, open_file, get_outdir, progress, check_indir, run_fasttree, run_iqtree, fix_iqtree

def main():

    # =============== #
    # 	PARAMETERS    #
    # =============== #

    args = docopt(__doc__,version='0.9.0')

    input_dir = check_indir(args['--input'])
    tree_results_file = args['--tree_results']
    mode = args['--mode']
    output_dir = get_outdir(args['--output'])
    threads = int(args['--threads'])

    mode_list = ['fasttree', 'iqtree']

    if mode not in mode_list:
        sys.exit("Mode " + mode + " is not valid. Please select from fasttree, iqtree")

    out_path = get_outdir(output_dir, add_dir="alt_topology")

    check_programs("FastTreeMP")
    check_programs("iqtree")
    print ("[STATUS]\t: Setting up")

    p = Pool(threads)
    fasttree_threads = 1
    iq_threads = 4

    # TOI, Gene, Exclude [1]
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

    for gene in tree:
        if str(type[gene]) == "HGT":

            selected_genes.append(gene)
            toi = ["EGP", "TOI"]
            toi.append(gene)
            filename = alignment[gene]

            toi_list = []
            free_list = []
            rest_list = []

            with open_file(filename) as handle:
                for record in SeqIO.parse(handle,"fasta"):
                    if any(n in record.id for n in toi):
                        toi_list.append(record.id)
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

            if mode == "fasttree":
                # Constaint file for FastTree
                total = len(toi_list) + len(free_list) + len(rest_list)
                con.write(" " + str(total) +" 1\n")
                for gene_l in toi_list:
                    con.write(gene_l + " 1\n")
                for gene_l in free_list:
                    con.write(gene_l + " -\n")
                for gene_l in rest_list:
                    con.write(gene_l + " 0\n")
                con.close()
                # FastTree params
                tree_path = os.path.join(out_path,gene + ".fasttree")
                fasttree_params = "-gamma -lg -constraints " + con_path + " " + alignment[gene] + " > " + tree_path
                t_list=[fasttree_params,fasttree_threads]
                fasttree_job_list.append(t_list)

                con_tree_list.append([tree_path,tree[gene],all_trees_path])
                #IqTree Alt Params
                iqmodel = "-m LG+G"
                iqtree_alt_params = iqmodel + " -n 0 -zb 10000 -au -pre " + tree_path + " -s " + alignment[gene] + " -z " + all_trees_path
                iqtree_alt_job_list.append(iqtree_alt_params)

            else:
                # Constaint file for IqTree
                toi_list = [w.replace('@', '_') for w in toi_list]
                rest_list = [w.replace('@', '_') for w in rest_list]
                con.write("((" + ', '.join(toi_list) + ")," + "" + ', '.join(rest_list) + ");" + "\n")
                con.close()
                #IqTree Params
                iqmodel = " -mset WAG,LG,JTT -AICc -mrate E,I,G,R "
                ufbootstrap = 1000
                tree_pre = os.path.join(out_path,gene)
                tree_path = tree_pre + ".treefile"
                iqtree_ml_params = "-quiet -nt " + str(iq_threads) + iqmodel + " -bb " + str(ufbootstrap) + " -g " + con_path + " -s " + alignment[gene] + " -pre " + tree_pre
                iqtree_ml_job_list.append(iqtree_ml_params)
                fix_iqtree_list.append(tree_path)

                con_tree_list.append([tree_path,tree[gene],all_trees_path])
                #IqTree Alt Params
                iqtree_alt_params = iqmodel + " -n 0 -zb 10000 -au -pre " + tree_path + " -s " + alignment[gene] + " -z " + all_trees_path
                iqtree_alt_job_list.append(iqtree_alt_params)

    if mode == "fasttree":
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

    for i, _ in enumerate(p.imap_unordered(run_iqtree,iqtree_alt_job_list),1):
        progress(i,1,len(iqtree_alt_job_list))

    alt_res_path = os.path.join(output_dir, "alt_topology_results.tsv")
    alt_res = open(alt_res_path,'w')
    alt_res.write("#gene"+ "\t" +"deltaL"+ "\t" + "bp-RELL" + "\t" +"p-KH"+ "\t" + "p-SH"+ "\t" + "c-ELW" + "\t" +"c-ELW"+ "\t" + "Significantly worse" + "\n")
    #Create final file
    for gene in selected_genes:
        if mode == "fasttree":
            p=subprocess.Popen("grep -P 'logL\s*del' -A 2 " + out_path + "/" + gene + "*fasttree.iqtree | tail -1", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        else:
            p=subprocess.Popen("grep -P 'logL\s*del' -A 2 " + out_path + "/" + gene + "*treefile.iqtree | tail -1", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        target_line = p.stdout.read().decode('utf').split()
        deltaL = target_line[2]
        ns = []
        pvalue = []
        for i in range(4,len(target_line),2):
            ns.append(target_line[i])
        significantly_worse = ns.count("-")
        for i in range(3,len(target_line),2):
            pvalue.append(target_line[i])
        alt_res.write(gene + "\t" + deltaL + "\t" + '\t'.join(pvalue) + "\t" + str(significantly_worse)+"\n")

    alt_res.close()



def concatenate_trees(con_tree_list):
    tree_path, tree_gene, all_trees_path = con_tree_list
    FNULL = open(os.devnull, 'w')
    subprocess.call("cat " + tree_path + " " + tree_gene + " | sed 's/@/_/g' | sed 's/|/_/g' > " + all_trees_path, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

if __name__ == '__main__':
    main()
