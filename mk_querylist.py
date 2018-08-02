import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from pandas import *
import pandas as pd
import textwrap
import time

def read_file(filename):
    try:
        usercols = ['desig_01', 'desig_02', 'ra', 'dec']
        sources = read_table(filename, header=0, names=usercols, comment='\\', sep='\s+',
                         usecols=[0, 1, 4, 5])
    except:
        raise IOError('The input file does not have the correct number of columns or they are in the incorrect order. Try again.')
    return sources

def mk_querylist(filename, output):

    # Read in the file
    sources = read_file(filename)

    # Initialize the output file
    if os.path.exists(output):
        print('{} already exists. Removing and Recreating'.format(output))
        os.remove(output)
    f = open(output, 'w')
    f.write("\ EQUINOX = 'J2000.0'\n")
    f.write("|         desig           |   ra      |  dec     |  major |\n")
    f.write("|         char            |   double  |  double  | double |\n")

    for i in range(len(sources.l)):
        f.write(' {} {}  {:10.6f} {:10.6f}      {:2.1f}\n'.format(sources.desig_01.iloc[i], sources.desig_02.iloc[i],
                                                                  sources.ra.iloc[i], sources.dec.iloc[i], 2.0))

    print('Done')

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
Description
------------------------------------

Locations
------------------------------------

Arguments
------------------------------------

Examples
------------------------------------
------------------------------------
''' ) )
    parser.add_argument("-f",
                        "--filename",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name and path of the archive or catalog to be processed",
                        default=None,
                        dest = "filename")
    parser.add_argument("-o",
                        "--output",
                        action='store',
                        nargs=1,
                        type=str,
                        help='The name of the output file',
                        default=None,
                        dest='output')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    print('args:', args)

    mk_querylist(args.filename[0], args.output[0])
