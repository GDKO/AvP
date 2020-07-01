#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Fix iqtree-generated trees
    Usage:
        fix_iqtree.py <tree>

"""
from docopt import docopt
import sys
import os
depot_path = os.path.join(sys.path[0], "../")
sys.path.insert(0,depot_path)


from depot.PetIO import fix_iqtree

args = docopt(__doc__)
fname = args['<tree>']

fix_iqtree(fname)
