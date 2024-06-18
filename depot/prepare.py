#!/usr/bin/env python3

"""
  Usage:
    avp prepare -a <FILE> -o <DIR> -f <FILE> -b <FILE> -x <FILE> -c <FILE> [-y <FILE>] [--cfg <FILE>]

  Options:
    -h, --help                  show this
    -a, --aifeatures <FILE>     alienness features file
    -o, --output <DIR>          creates a directory for all output files
    -f, --fastafile <FILE>      fasta file of proteome
    -b, --hitsfile <FILE>      hits file in -outfmt 6
    -x, --tax_groups <FILE>     file specifying taxonomic groups (ingroup, exclude)
    -c, --config_file <FILE>    config file
    -y, --ortho_groups <FILE>   Group genes based on an orthology groups file
    --cfg <FILE>                Specify a taxonomy config file (default file at depot/taxonomy.yaml)
"""

import sys
import os
import subprocess
import csv
import re
import numpy as np
import networkx as nx
import yaml

from ete3 import NCBITaxa
from Bio import SeqIO
from docopt import docopt
from multiprocessing import Pool
from operator import itemgetter

from depot.PetIO import check_programs, get_outdir, progress, open_file

def main():

    # =============== #
    # 	PARAMETERS    #
    # =============== #

    args = docopt(__doc__)

    ai_features = args['--aifeatures']
    output_dir = get_outdir(args['--output'])
    fasta_inputfile = args['--fastafile']
    hits_inputfile = args['--hitsfile']
    groups_yaml = args['--tax_groups']
    config_yaml = args['--config_file']
    ortho_groups = args['--ortho_groups']
    cfg_file = args['--cfg']

    # =============== #
    # 	MAIN          #
    # =============== #

    """
    0. Setting up
    """

    print ("[+] Setting up")

    # Check if programs in path
    check_programs("blastdbcmd", "mafft")

    # Create folders
    fasta_folder = os.path.join(output_dir,"fastagroups")
    get_outdir(fasta_folder)

    mafft_folder = os.path.join(output_dir,"mafftgroups")
    get_outdir(mafft_folder)

    tmp_folder = os.path.join(output_dir,"tmp")
    get_outdir(tmp_folder)

    # Load proteome in memory
    record_dict = SeqIO.index(fasta_inputfile, "fasta")

    # Create taxonomic groups
    orgtag = "@StudiedOrganism"

    stream = open(groups_yaml, 'r')
    ingroup_egp = yaml.safe_load(stream)
    stream.close()

    if not cfg_file:
        cfg_file = os.path.join(sys.path[0],"depot","taxonomy.yaml")
    stream = open(cfg_file, 'r')
    config_groups = yaml.safe_load(stream)
    stream.close()

    stream = open(config_yaml, 'r')
    config_opts = yaml.safe_load(stream)
    stream.close()

    data_type = config_opts["data_type"]
    threads = config_opts["max_threads"]
    trim = config_opts["trimal"]
    ai_cutoff = config_opts["ai_cutoff"]
    ahs_cutoff = config_opts["ahs_cutoff"]
    pct_cutoff = config_opts["outg_pct_cutoff"]
    selection = config_opts["selection"]
    percent_identity = config_opts["percent_identity"]
    cutoffextend = config_opts["cutoffextend"]
    min_num_hits = config_opts["min_num_hits"]
    percentage_similar_hits  = config_opts["percentage_similar_hits"]
    mode = config_opts["mode"]
    mafft_options = config_opts["mafft_options"]
    trimal_options = config_opts["trimal_options"]

    # turn selection into condition
    selection = selection.replace("ai","float(row[i_ai])>"+str(ai_cutoff))
    selection = selection.replace("ahs","float(row[i_ahs])>"+str(ahs_cutoff))
    selection = selection.replace("outg_pct","float(row[i_pct])>"+str(pct_cutoff))

    forbidden_chars = ["|", "@", ":"]

    if trim:
        check_programs("trimal")
        trim_folder = os.path.join(output_dir,"trim")
        get_outdir(trim_folder)
    else:
        trim_folder = ""

    dbtype = ""

    if data_type == "AA":
        dbtype = "prot"
    elif data_type == "DNA":
        dbtype = "nucl"
    else:
        sys.exit("[x] data_type should be either AA or DNA")

    #Setting up NCBI Taxonomy
    ncbi = NCBITaxa()

    """
    1. Select HGT
    """

    query_dict_set = {}
    queries_info = {}

    with open(ai_features, 'r', encoding='utf8') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')

        for row in reader:
            # Get index col
            i_query = row.index('query name')
            i_donor = row.index('donor')
            i_ingroup = row.index('ingroup')
            i_ai = row.index('AI')
            i_hgt = row.index('HGTindex')
            i_nbhits = row.index('query hits number')
            i_ahs = row.index('AHS')
            i_pct = row.index('outg_pct')
            break

        for row in reader:
            L_donor = row[i_donor].rstrip('\n').rsplit(':',4)
            L_ingroup = row[i_ingroup].rstrip('\n').rsplit(':',4)
            if (row[i_donor] != '::::'): #Skipping hits with only Ingroup
                if (float(row[i_nbhits])>= min_num_hits and float(L_donor[2]) <= percent_identity):
                    if eval(selection):
                        donor_pos = int(L_donor[1])
                        if(row[i_ingroup] == '::::'):
                            ingroup_pos = 0
                        else:
                            ingroup_pos = int(L_ingroup[1])
                        #Select at least 50 hits
                        last_pos = min(max(max(ingroup_pos,donor_pos) + cutoffextend, 50), int(row[i_nbhits]))
                        queries_info[row[i_query]] = {'pos':last_pos}
                        query_dict_set[row[i_query]] = set()

    print ("[!] Selected " + str(len(query_dict_set)) + " HGT candidates")

    """
    2. Parse hits file
    """

    print ("[+] Parsing hits file and grouping similar queries")

    extract_hit_id_set = set()

    with open_file(hits_inputfile) as fhr_bl:
        for line in fhr_bl:
            if('#' not in line):
                L_hitqline = line.rstrip('\n').split('\t')
                query_id = L_hitqline[0]
                if any(char in query_id for char in forbidden_chars):
                    sys.exit("[x] Query Ids should not contain | or @ or : characters")
                if len(query_id) > 255:
                    sys.exit("[x] Query names are too long")
                if query_id in queries_info.keys(): # Queries that pass the initial selection
                    if(len(query_dict_set[query_id]) <= queries_info[query_id]["pos"]):
                        query_hit_id = L_hitqline[1]
                        if "@" in query_id:
                            sys.exit("@ symbol is not allowed: " + query_id)
                        if "@" in query_hit_id:
                            sys.exit("@ symbol is not allowed: " + query_hit_id)
                        extract_hit_id_set.add(query_hit_id)
                        query_dict_set[query_id].add(query_hit_id) # GK

    # Group hits
    G = nx.Graph()

    if ortho_groups:
        num_groups = 0
        with open_file(ortho_groups) as fhr_og:
            for line in fhr_og:
                num_groups += 1
                members = line.split()
                for i in range(1,len(members),1):
                    G.add_node(members[i])
                    if i>1:
                        G.add_edge(members[i], members[i-1])
        print ("[!] Found " + str(num_groups) + " groups")

    else:
        for protein_id, hits in query_dict_set.items():
            G.add_node(protein_id)
            for protein_id_other, hitsc in query_dict_set.items():
                if protein_id != protein_id_other:
                    u = len(set.intersection(hits, hitsc))
                    m = min(len(hits),len(hitsc))
                    if (u/m) >= percentage_similar_hits:
                        G.add_edge(protein_id, protein_id_other)
        print ("[!] Formed " + str(len(list(nx.connected_components(G)))) + " groups")

    """
    3. Extract hits
    """

    print ("[+] Extracting hits from DB")

    extract_id_path = os.path.join(tmp_folder,"extract_id.txt")
    fhw_extract_id = open(extract_id_path, 'w')
    fhw_extract_id.write('\n'.join(extract_hit_id_set)+'\n')
    fhw_extract_id.close()

    setblastfa_path = os.path.join(tmp_folder,"setblast.fa")
    fhw_setblastfa = open(setblastfa_path, 'w')
    setblastlog_path = os.path.join(tmp_folder,"setblast.log")

    if mode == "blast":
        blastdbcmd_command = 'blastdbcmd -db '+ config_opts["blast_db_path"] + ' -dbtype ' + dbtype + ' -entry_batch ' +  extract_id_path + ' -outfmt ">%a@%T\n%s" -logfile ' + setblastlog_path + ' -out ' + setblastfa_path
        subprocess.call(blastdbcmd_command, shell= True)
    else: # GK This is specific to SwissProt for now, have to test for UniProt in the future
        if mode == "sp":
            db_re = re.compile("OX=[0-9]*")
        elif mode == "ur90":
            db_re = re.compile("TaxID=[0-9]*")
        else:
            sys.exit(mode + " is not a valid mode")
        with open_file(config_opts["fasta_path"]) as handle:
            for record in SeqIO.parse(handle,"fasta"):
                if record.id in extract_hit_id_set:
                    ox, taxid = db_re.search(record.description).group().split("=")
                    fhw_setblastfa.write(">" + record.id + "@" + taxid + "\n")
                    seq = str(record.seq)
                    fhw_setblastfa.write(seq + "\n")

    fhw_setblastfa.close()

    # Load hits to memory
    hits_dict = {}
    record_to_taxid = {}

    with open_file(setblastfa_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id,taxid = record.id.rstrip('\n').split('@')
            record_to_taxid[id] = taxid
            hits_dict[id] = str(record.seq)

    """
    4. Write fasta
    """
    print ("[+] Writing fasta files")

    taxonomy_nexus_path = os.path.join(tmp_folder,"taxonomy_nexus.txt")
    taxonomy_nexus = open(taxonomy_nexus_path,'w')
    groups_tsv_path = os.path.join(output_dir,"groups.tsv")
    groups_tsv = open(groups_tsv_path,'w')

    group_id = 1
    final_number_of_candidates = 0
    final_number_of_groups = 0

    group_dict = {}
    number_of_lost_taxids = 0
    number_of_lost_records = 0

    for subgraph in nx.connected_components(G):

        num_of_seqs_in_group = 0
        grouped_hits = []
        queries_to_fasta = []
        queries_to_group = []

        for q in subgraph:
            queries_to_fasta.append(q)
            queries_to_group.append(q)
            grouped_hits.extend(query_dict_set[q])

        grouped_hits = set(grouped_hits)

        if len(grouped_hits) >= min_num_hits:
            group_name_file = "gp" + str(group_id) + '.fa'
            gp_pathname = os.path.join(fasta_folder,group_name_file)
            fw_gp = open(gp_pathname,'w')
            final_number_of_groups += 1

            for q in queries_to_fasta:
                final_number_of_candidates += 1
                fw_gp.write('>' + q + orgtag + '\n' + str(record_dict[q].seq) + '\n')
                num_of_seqs_in_group += 1

            for record_id in grouped_hits:
                if record_id in hits_dict:
                    taxid = record_to_taxid[record_id]
                    try:
                        ncbi.get_lineage(taxid)
                        taxid_found = True
                    except:
                        taxid_found = False
                    if not taxid or not taxid_found:
                        # actually print a file containing lost taxids
                        number_of_lost_taxids += 1
                        selectname = "Unknown"
                    else:
                        lnode = set(ncbi.get_lineage(taxid))
                        lname = ncbi.translate_to_names(lnode) # This does not output them in order
                        #print(ncbi.get_rank(lnode)) Maybe try this with {1: 'no rank', 2: 'superkingdom'}
                        taxonomy_nexus.write(str(record_id)+"\t"+str(lname)+"\n")

                        egp_hit = list(lnode.intersection(set(ingroup_egp["EGP"].keys())))
                        ingroup_hit = list(lnode.intersection(set(ingroup_egp["Ingroup"].keys())))
                        cfg_hit = list(lnode.intersection(set(config_groups["Other"].keys())))
                        kdom_hit = list(lnode.intersection(set(config_groups["Kingdom"].keys())))

                        if egp_hit:
                            selectname = "EGP-" + ingroup_egp["EGP"][egp_hit[0]]
                        elif ingroup_hit:
                            selectname = "Ingroup-" + ingroup_egp["Ingroup"][ingroup_hit[0]]
                        elif cfg_hit:
                            selectname = config_groups["Other"][cfg_hit[0]]
                        elif kdom_hit:
                            selectname = config_groups["Kingdom"][kdom_hit[0]]
                        else:
                            selectname = "Unknown"

                    fw_gp.write(">" + record_id + "@" + selectname + "\n")
                    fw_gp.write(hits_dict[record_id] + "\n")
                    num_of_seqs_in_group += 1
                else:
                    # actually print a file containing lost gids
                    number_of_lost_records += 1

            groups_tsv.write(group_name_file + '\t' + str(num_of_seqs_in_group) + '\t' + '\t'.join(queries_to_group) + '\n')
            group_dict[group_name_file] = num_of_seqs_in_group
            fw_gp.close()
            group_id += 1

    groups_tsv.close()

    print ("[!] Skipped " + str(number_of_lost_records) + " hits and " + str(number_of_lost_taxids) + " taxids.")

    """
    5. Align fasta
    """

    print ("[+] Aligning fasta files")

    jobs = threads
    p = Pool(jobs)

    job_list = []

    for group_name, value in sorted(group_dict.items(), key = itemgetter(1), reverse = True):
        g_list = [group_name,fasta_folder,mafft_folder,trim,trim_folder, mafft_options, trimal_options]
        job_list.append(g_list)

    i = 0
    for i, _ in enumerate(p.imap_unordered(run_mafft, job_list),1):
        progress(i,1,len(job_list))

    print ("[!] Finished with " + str(final_number_of_candidates) + " HGT candidates in " + str(final_number_of_groups) + " groups")


def run_mafft(g_list):
    group_name = g_list[0]
    fasta_folder = g_list[1]
    mafft_folder = g_list[2]
    trim = g_list[3]
    trim_folder = g_list[4]
    mafft_options = g_list[5]
    trimal_options = g_list[6]
    group_fasta_file_path = os.path.join(fasta_folder,group_name)
    group_align_file_path = os.path.join(mafft_folder,group_name)
    FNULL = open(os.devnull, 'w')
    subprocess.call("mafft " + mafft_options + " --thread 1 " + group_fasta_file_path + " > " + group_align_file_path, shell=True, stderr=FNULL)

    if trim:
        group_trim_file_path = os.path.join(trim_folder,group_name)
        subprocess.call("trimal " + trimal_options + " -in " + group_align_file_path + " -out " + group_trim_file_path, shell=True, stderr=FNULL)

if __name__ == '__main__':
    main()
