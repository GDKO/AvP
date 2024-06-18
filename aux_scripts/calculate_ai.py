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
    outout_file_path = input_file + "_ai.out"
    outout_file = open(outout_file_path,'w')

    stream = open(groups_yaml, 'r')
    ingroup_egp = yaml.safe_load(stream)
    stream.close()

    ingroup_taxid=set(ingroup_egp["Ingroup"].keys())
    egp_taxid=set(ingroup_egp["EGP"].keys())
    egp_taxid.add(2787823)
    egp_taxid.add(2787854)

    forbidden_chars = ["|", "@", ":"]

    #Setting up NCBI Taxonomy
    ncbi = NCBITaxa()
    number_of_lost_taxids = 0

    list_genes = []
    best_hit_ingroup = {}
    best_hit_donor = {}
    num_hits = {}
    sum_ingroup_bitscore = {} #GK_ahs
    sum_donor_bitscore = {} #GK_ahs
    set_taxid_ingroup = {} #GK_outg_pct
    set_taxid_donor = {} #GK_outg_pct


    outout_file.write("query name\tdonor\tingroup\tAI\tHGTindex\tquery hits number\tAHS\toutg_pct\n")

    skip = 0
    with open_file(input_file) as fhr_bl:
        for line in fhr_bl:
            elements = line.split()
            if len(elements)<13:
                skip += 1
            else:
                gene = elements[0]
                if any(char in gene for char in forbidden_chars):
                    sys.exit("[x] Query Ids should not contain | or @ or : characters")
                if len(gene) > 255:
                    sys.exit("[x] Query names are too long")
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
                    sum_ingroup_bitscore[gene] = [] #GK_ahs
                    sum_donor_bitscore[gene] = [] #GK_ahs
                    set_taxid_ingroup[gene] = set() #GK_outg_pct
                    set_taxid_donor[gene] = set() #GK_outg_pct
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
                    elif lnode.intersection(ingroup_taxid):
                        sum_ingroup_bitscore[gene].append(float(bitscore)) #GK_ahs
                        set_taxid_ingroup[gene].add(taxid) #GK_outg_pct
                        if gene not in best_hit_ingroup.keys():
                            best_hit_ingroup[gene] = {}
                            best_hit_ingroup[gene]["hit"] = hit
                            best_hit_ingroup[gene]["pos"] = str(num_hits[gene])
                            best_hit_ingroup[gene]["iden"] = iden
                            best_hit_ingroup[gene]["evalue"] = str(evalue)
                            best_hit_ingroup[gene]["bitscore"] = bitscore
                    else:
                        sum_donor_bitscore[gene].append(float(bitscore)) #GK_ahs
                        set_taxid_donor[gene].add(taxid) #GK_outg_pct
                        if gene not in best_hit_donor.keys():
                            best_hit_donor[gene] = {}
                            best_hit_donor[gene]["hit"] = hit
                            best_hit_donor[gene]["pos"] = str(num_hits[gene])
                            best_hit_donor[gene]["iden"] = iden
                            best_hit_donor[gene]["evalue"] = str(evalue)
                            best_hit_donor[gene]["bitscore"] = bitscore

    for gene in set(list_genes):

        if gene not in best_hit_ingroup.keys():
            ingroup_str = "::::"
            ingroup_evalue = 1
            ingroup_bitscore = 0
        else:
            ingroup_str = best_hit_ingroup[gene]["hit"] + ":" + best_hit_ingroup[gene]["pos"] + ":" + best_hit_ingroup[gene]["iden"] + ":" + best_hit_ingroup[gene]["evalue"] + ":" + best_hit_ingroup[gene]["bitscore"]
            ingroup_evalue = float(best_hit_ingroup[gene]["evalue"])
            ingroup_bitscore = float(best_hit_ingroup[gene]["bitscore"])

        if gene not in best_hit_donor.keys():
            donor_str = "::::"
            donor_evalue = 1
            donor_bitscore = 0
        else:
            donor_str = best_hit_donor[gene]["hit"] + ":" + best_hit_donor[gene]["pos"] + ":" + best_hit_donor[gene]["iden"] + ":" + best_hit_donor[gene]["evalue"] + ":" + best_hit_donor[gene]["bitscore"]
            donor_evalue = float(best_hit_donor[gene]["evalue"])
            donor_bitscore = float(best_hit_donor[gene]["bitscore"])

        ai = calculate_ai(ingroup_evalue,donor_evalue)
        hgt_score = donor_bitscore - ingroup_bitscore
        if len(set_taxid_donor[gene]) == 0:
            outg_pct = 0 #GK_outg_pct
        else:
            outg_pct = round(100 * len(set_taxid_donor[gene]) / (len(set_taxid_donor[gene]) + len(set_taxid_ingroup[gene]))) #GK_outg_pct
        sum_ingroup_bitscore[gene].sort(reverse=True)
        sum_donor_bitscore[gene].sort(reverse=True)

        norm_bit_ingroup_gene = 0
        norm_bit_donor_gene = 0
        for bit in sum_ingroup_bitscore[gene]:
            norm_bit = calculate_norm_bitscore(bit,sum_ingroup_bitscore[gene][0],a)
            norm_bit_ingroup_gene += norm_bit
        for bit in sum_donor_bitscore[gene]:
            norm_bit = calculate_norm_bitscore(bit,sum_donor_bitscore[gene][0],a)
            norm_bit_donor_gene += norm_bit
        Dnorm = norm_bit_donor_gene - norm_bit_ingroup_gene
        #Dchs = sum(sum_donor_bitscore[gene]) - sum(sum_ingroup_bitscore[gene])
        outout_file.write(gene + "\t" + donor_str + "\t" + ingroup_str + "\t" + str(ai) + "\t" + str(hgt_score) + "\t" + str(num_hits[gene])+"\t")
        outout_file.write(str(Dnorm) + "\t" + str(outg_pct) + "\n")

    outout_file.close()

    print("[!] Skipped " + str(skip) + " hits") 

def calculate_ai(ingroup_evalue,donor_evalue):
    offset = 1e-200
    ai = math.log(ingroup_evalue + offset) - math.log(donor_evalue + offset);
    return ai

def calculate_norm_bitscore(bitscore,hbitscore,a):
    norm_bitscore = bitscore*math.exp(-a*((hbitscore-bitscore)/hbitscore))
    return norm_bitscore

if __name__ == '__main__':
    main()
