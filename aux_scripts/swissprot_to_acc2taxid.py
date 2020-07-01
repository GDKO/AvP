#!/usr/bin/env python3

"""swissprot_to_acc2taxid.py
  Usage:
    swissprot_to_acc2taxid.py -i <FILE>

  Options:
    -h, --help                    show this
    -i, --input <FILE>             folder with aligned fasta files
    --version                     print version
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
    args = docopt(__doc__,version='0.9.0')
    input_file = args['--input']


    print("acc\tver\ttaxid\n")

    swissprot_re = re.compile("OX=\d*")
    swissprot_re_sv = re.compile("SV=\d*")

    with open_file(input_file) as handle:
        for record in SeqIO.parse(handle,"fasta"):
            ver = ""
            id = ""
            taxid = ""
            tr = ""
            tr, id, tr = record.name.split("|")
            print(id)
            tr, taxid = swissprot_re.search(record.description).group().split("=")
            tr, ver = swissprot_re_sv.search(record.description).group().split("=")



if __name__ == '__main__':
    main()
