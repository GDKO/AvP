#!/usr/bin/env python3

"""calculate_ai.py
    Usage:
        calculate_ai.py -i <FILE> -x <FILE>

    Options:
        -i, --input <FILE>          outfmt format 6 file in the following format std + taxid
        -x, --tax_groups <FILE>     file specifying taxonomic groups (ingroup, exclude)
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
    outout_file_path = input_file + "_ai.out"
    outout_file = open(outout_file_path,'w')

    stream = open(groups_yaml, 'r')
    toi_egp = yaml.safe_load(stream)
    stream.close()

    toi_taxid=set(toi_egp["TOI"].keys())
    egp_taxid=set(toi_egp["EGP"].keys())

    #Setting up NCBI Taxonomy
    ncbi = NCBITaxa()
    number_of_lost_taxids = 0

    list_genes = []
    best_hit_toi = {}
    best_hit_ntoi = {}
    num_hits = {}

    outout_file.write("query name\tdonor\trecipient\tAI\tHGTindex\tquery hits number\n")

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
                    if gene not in best_hit_toi.keys():
                        best_hit_toi[gene] = {}
                        best_hit_toi[gene]["hit"] = hit
                        best_hit_toi[gene]["pos"] = str(num_hits[gene])
                        best_hit_toi[gene]["iden"] = iden
                        best_hit_toi[gene]["evalue"] = str(evalue)
                        best_hit_toi[gene]["bitscore"] = bitscore
                else:
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
        outout_file.write(gene + "\t" + ntoi_str + "\t" + toi_str + "\t" + str(ai) + "\t" + str(hgt_score) + "\t" + str(num_hits[gene])+"\n")

    outout_file.close()

def calculate_ai(toi_evalue,ntoi_evalue):
    offset = 1e-200
    ai = math.log(toi_evalue + offset) - math.log(ntoi_evalue + offset);
    return ai


if __name__ == '__main__':
    main()
