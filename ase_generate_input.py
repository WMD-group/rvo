#! /usr/bin/env python

"""
Handy little script for generating RVO input. Requires a set of output files compatible with the Atomic Simulation Environment.

"""

import ase.io
import argparse

def main(paths):
    for path in paths:
        A = ase.io.read(path)
        cell = A.get_cell()
        print A.get_total_energy(),
        for x in cell.flatten():
            print ',{0}'.format(x),
        print ''

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description='Generate 10-column input file for RVO from given quantum chemistry calculation files')
    parser.add_argument('paths', nargs='+')
    args = parser.parse_args()
    main(args.paths)
