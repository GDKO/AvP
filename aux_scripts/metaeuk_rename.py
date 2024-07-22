#!/usr/bin/env python3

"""metaeuk_rename.py
    Usage:
        metaeuk_rename.py -f <FILE> -g <FILE> [-p <STR>]

    Options:
        -f, --fasta <FILE>      protein fasta file
        -g, --gff <FILE>        gff file
        -p, --prefix <STR>      prefix for headers [default: protein]
        -h, --help              show this
"""

import sys
import os
depot_path = os.path.join(sys.path[0], "../")
sys.path.insert(0,depot_path)

from docopt import docopt
from depot.PetIO import open_file
from Bio import SeqIO

def main():
    args = docopt(__doc__)

    fasta_file = args['--fasta']
    gff_file = args['--gff']
    prefix = args['--prefix']

    output_fasta_file_path = fasta_file + "_renamed.faa"
    output_gff_file_path = gff_file + "_renamed.bed"

    header_table = {}

    i = 1

    with open_file(fasta_file) as handle, open(output_fasta_file_path, "w", encoding = "UTF-8") as fo:
        for record in SeqIO.parse(handle, "fasta"):
            header_lst = record.id.split("|")
            new_header = prefix + "_" + str(i)
            old_header = header_lst[0] + "|" + header_lst[1] + "|" + header_lst[2] + "|" + header_lst[6]
            header_table[old_header] = new_header
            fo.write(">" + new_header + "\n" + str(record.seq) + "\n")
            i += 1 
    
    with open_file(gff_file) as handle, open(output_gff_file_path, "w", encoding = "UTF-8") as fo:
        for line in handle:
            gff_lst = line.rstrip().split()
            old_header = line.rstrip().split("=")[-1]
            if gff_lst[2] == "mRNA":
                fo.write(gff_lst[0] + "\t" + gff_lst[3] + "\t" + gff_lst[4] + "\t" + header_table[old_header] + "\n")



if __name__ == '__main__':
    main()
