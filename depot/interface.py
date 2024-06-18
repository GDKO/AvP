#!/usr/bin/env python3

"""
  Usage:
    avp [<module>] [<args>...] [-h] [-v]

  Modules:
    prepare         extract the candidates, prepare the groups, and perform MSA
    detect          perform tree inference and detect HGT candidates
    classify        classify genes in specific taxonomic categories

    evaluate        evaluate alternative topology

  Options
    -h, --help      show this
    -v, --version   show version number

    See 'avp <command> --help' for more information on a specific command.

"""

import sys
from docopt import docopt

def main():
    args = docopt(__doc__,version='1.0.5', options_first=True)
    if args['<module>'] == 'prepare':
        import depot.prepare as prepare
        prepare.main()
    elif args['<module>'] == 'detect':
        import depot.detect as detect
        detect.main()
    elif args['<module>'] == 'classify':
        import depot.classify as classify
        classify.main()
    elif args['<module>'] == 'evaluate':
        import depot.evaluate as evaluate
        evaluate.main()
    else:
        sys.exit("%r is not an avp module. See 'avp -h'." % args['<module>'])
