#!/usr/bin/env python3

"""
  Usage:
    acc2taxid.py -i <FILE> -m <STR>

  Options:
    -h, --help                    show this
    -i, --input <FILE>            folder with aligned fasta files
    -m, --mode <STR>              mode, choose between swissprot or uniref
"""

import sys
import os
depot_path = os.path.join(sys.path[0], "../")
sys.path.insert(0,depot_path)

import re
from Bio import SeqIO
from docopt import docopt
from depot.PetIO import open_file

def main():
    args = docopt(__doc__)
    input_file = args['--input']
    mode = args['--mode']


    if mode == "swissprot":
        db_re = re.compile("OX=[0-9]*")
    elif mode == "uniref":
        db_re = re.compile("TaxID=[0-9]*")
    else:
        sys.exit(mode + " is not a valid mode")

    print("accession.version\ttaxid")

    with open_file(input_file) as handle:
        for record in SeqIO.parse(handle,"fasta"):
            if mode == "swissprot":
                tr, id, tr = record.name.split("|")
            elif mode == "uniref":
                tr, id = record.name.split("|")
            tr, taxid = db_re.search(record.description).group().split("=")
            if not taxid:
                taxid = 1
            print(id + "\t" + str(taxid))

if __name__ == '__main__':
    main()
