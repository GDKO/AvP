#!/usr/bin/env python3

"""calculate_ai.py
    Usage:
        calculate_ai.py -i <FILE> -x <FILE> [-a <INT>]

    Options:
        -i, --input <FILE>          outfmt format 6 file in the following format std + taxid
        -x, --tax_groups <FILE>     file specifying taxonomic groups (ingroup, exclude)
        -a, --alpha <INT>           alpha parameter for the normalised score [default: 10]
        -h, --help                  show this
"""

import sys
import os
depot_path = os.path.join(sys.path[0], "../")
sys.path.insert(0,depot_path)

import yaml
import math
from ete3 import NCBITaxa
from docopt import docopt
from depot.PetIO import open_file

def main():
    args = docopt(__doc__)

    input_file = args['--input']
    groups_yaml = args['--tax_groups']
    a = int(args['--alpha'])
    chs_threshold = 0.8 #GK_chs
    outout_file_path = input_file + "_ai.out"
    outout_file = open(outout_file_path,'w')

    stream = open(groups_yaml, 'r')
    toi_egp = yaml.safe_load(stream)
    stream.close()

    toi_taxid=set(toi_egp["TOI"].keys())
    egp_taxid=set(toi_egp["EGP"].keys())
    egp_taxid.add(2787823)
    egp_taxid.add(2787854)

    #Setting up NCBI Taxonomy
    ncbi = NCBITaxa()
    number_of_lost_taxids = 0

    list_genes = []
    best_hit_toi = {}
    best_hit_ntoi = {}
    num_hits = {}
    sum_toi_bitscore = {} #GK_chs
    sum_ntoi_bitscore = {} #GK_chs
    num_ntoi = {} #GK_chs
    num_toi = {} #GK_chs

    outout_file.write("query name\tdonor\trecipient\tAI\tHGTindex\tquery hits number\tAHS\n")

    with open_file(input_file) as fhr_bl:
        for line in fhr_bl:
            elements = line.split()
            gene = elements[0]
            list_genes.append(gene)
            hit = elements[1]
            iden = elements[2]
            evalue = elements[10]
            bitscore = elements[11]
            taxid = elements[-1]
            if ";" in taxid:
                first_id = taxid.split(";")[0] #get the first id from taxid if multiple
                taxid = first_id
            if gene not in num_hits.keys():
                num_hits[gene] = 0
                sum_toi_bitscore[gene] = [] #GK_chs
                sum_ntoi_bitscore[gene] = [] #GK_chs
                num_ntoi[gene] = 0 #GK_chs
                num_toi[gene] = 0 #GK_chs
            num_hits[gene] += 1
            try:
                ncbi.get_lineage(taxid)
            except:
                # actually print a file containing lost taxids
                number_of_lost_taxids += 1
            else:
                lnode = set(ncbi.get_lineage(taxid))
                if lnode.intersection(egp_taxid):
                    egp = 1
                elif lnode.intersection(toi_taxid):
                    sum_toi_bitscore[gene].append(float(bitscore)) #GK_chs
                    num_toi[gene] += 1 #GK_chs
                    if gene not in best_hit_toi.keys():
                        best_hit_toi[gene] = {}
                        best_hit_toi[gene]["hit"] = hit
                        best_hit_toi[gene]["pos"] = str(num_hits[gene])
                        best_hit_toi[gene]["iden"] = iden
                        best_hit_toi[gene]["evalue"] = str(evalue)
                        best_hit_toi[gene]["bitscore"] = bitscore
                else:
                    sum_ntoi_bitscore[gene].append(float(bitscore)) #GK_chs
                    num_ntoi[gene] += 1 #GK_chs
                    if gene not in best_hit_ntoi.keys():
                        best_hit_ntoi[gene] = {}
                        best_hit_ntoi[gene]["hit"] = hit
                        best_hit_ntoi[gene]["pos"] = str(num_hits[gene])
                        best_hit_ntoi[gene]["iden"] = iden
                        best_hit_ntoi[gene]["evalue"] = str(evalue)
                        best_hit_ntoi[gene]["bitscore"] = bitscore

    for gene in set(list_genes):

        if gene not in best_hit_toi.keys():
            toi_str = "::::"
            toi_evalue = 1
            toi_bitscore = 0
        else:
            toi_str = best_hit_toi[gene]["hit"] + ":" + best_hit_toi[gene]["pos"] + ":" + best_hit_toi[gene]["iden"] + ":" + best_hit_toi[gene]["evalue"] + ":" + best_hit_toi[gene]["bitscore"]
            toi_evalue = float(best_hit_toi[gene]["evalue"])
            toi_bitscore = float(best_hit_toi[gene]["bitscore"])

        if gene not in best_hit_ntoi.keys():
            ntoi_str = "::::"
            ntoi_evalue = 1
            ntoi_bitscore = 0
        else:
            ntoi_str = best_hit_ntoi[gene]["hit"] + ":" + best_hit_ntoi[gene]["pos"] + ":" + best_hit_ntoi[gene]["iden"] + ":" + best_hit_ntoi[gene]["evalue"] + ":" + best_hit_ntoi[gene]["bitscore"]
            ntoi_evalue = float(best_hit_ntoi[gene]["evalue"])
            ntoi_bitscore = float(best_hit_ntoi[gene]["bitscore"])

        ai = calculate_ai(toi_evalue,ntoi_evalue)
        hgt_score = ntoi_bitscore - toi_bitscore
        if num_ntoi[gene] == 0:
            per_ntoi = 0
        else:
            per_ntoi = num_ntoi[gene]/(num_ntoi[gene] + num_toi[gene]) #GK_chs
        sum_toi_bitscore[gene].sort(reverse=True)
        sum_ntoi_bitscore[gene].sort(reverse=True)

        if (sum(sum_ntoi_bitscore[gene]) > sum(sum_toi_bitscore[gene]) and per_ntoi>chs_threshold): #GK_chs
            chs = "Y" #GK_chs
        else: #GK_chs
            chs = "N" #GK_chs
        norm_bit_toi_gene = 0
        norm_bit_ntoi_gene = 0
        for bit in sum_toi_bitscore[gene]:
            norm_bit = calculate_norm_bitscore(bit,sum_toi_bitscore[gene][0],a)
            norm_bit_toi_gene += norm_bit
        for bit in sum_ntoi_bitscore[gene]:
            norm_bit = calculate_norm_bitscore(bit,sum_ntoi_bitscore[gene][0],a)
            norm_bit_ntoi_gene += norm_bit
        Dnorm = norm_bit_ntoi_gene - norm_bit_toi_gene
        Dchs = sum(sum_ntoi_bitscore[gene]) - sum(sum_toi_bitscore[gene])
        outout_file.write(gene + "\t" + ntoi_str + "\t" + toi_str + "\t" + str(ai) + "\t" + str(hgt_score) + "\t" + str(num_hits[gene])+"\t")
        #outout_file.write(str(Dchs) + "\t" + str(per_ntoi) + "\t" + chs + "\t" + str(Dnorm) + "\n") #GK_chs
        outout_file.write(str(Dnorm) + "\n")

    outout_file.close()

def calculate_ai(toi_evalue,ntoi_evalue):
    offset = 1e-200
    ai = math.log(toi_evalue + offset) - math.log(ntoi_evalue + offset);
    return ai

def calculate_norm_bitscore(bitscore,hbitscore,a):
    norm_bitscore = bitscore*math.exp(-a*((hbitscore-bitscore)/hbitscore))
    return norm_bitscore

if __name__ == '__main__':
    main()
