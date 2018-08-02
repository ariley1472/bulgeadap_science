#! /usr/bin/env python

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
import argparse
import glob
import numpy as np
from math import sqrt
import os
from pandas import *
import pandas as pd
import scipy.spatial as spatial
import textwrap


def convert_to_arcsec(long, lat):
    try:
        l1 = Angle(long, u.deg)
        b1 = Angle(lat, u.deg)
    except:
        print(long, lat)
        raise KeyboardInterrupt

    return [l1.arcsec, b1.arcsec]


def read_file(filename):
    # Use Pandas to read in the file
    try:
        usercols = ['desig_01', 'desig_02', 'l', 'b', 'ra', 'dec', 'magj', 'dmagj', 'magh', 'dmagh', 'magk',
                    'dmagk', 'mag3', 'dmag3', 'mag4', 'dmag4', 'mag5', 'dmag5', 'mag8', 'dmag8']
        sources = read_table(filename, header=3, names=usercols, comment='\\', sep='\s+',
                             usecols=[0, 1, 4, 5, 8, 9, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26])

    except:
        usercols = ['desig_01', 'desig_02']
        sources = read_table(filename, header=1, names=usercols, comment='\\', sep='\s+',
                             usecols=[0, 1])


    print(sources)
    raise KeyboardInterrupt
    return sources


def find_matches(epoch1_coords, deep_coords):
    allepoch_matches_ind = []
    ind = []
    count = 0

    epoch1_cat = SkyCoord(deep_coords[0], deep_coords[1], unit='deg', frame='galactic')
    for i in range(len(deep_coords[0])):
        print('i = ', i)
        target = SkyCoord([deep_coords[0][i], ], [deep_coords[1][i], ], unit='deg', frame='galactic')
        idxsearcharound, epoch1_match, sep2d, dist3d = epoch1_cat.search_around_sky(target, 2 * u.arcsec)

        if len(epoch1_match) == 0:  # No matches
            print('No Match Here!')
        elif len(epoch1_match) == 1:  # One match
            print('We Have A Match!')
            count += 1
            deep_matches_ind.append(i)

    print('length:', len(deep_matches_ind))
    raise KeyboardInterrupt
    return deep_matches_ind


def catalog_diff(epoch1_file, output):
    # read in all the necessary files:
    e1_sources = read_file(epoch1_file)
    e1_l = []
    e1_b = []
    for i in range(len(e1_sources)):
        #print(e1_sources.desig_02[i])
        e1_l.append(float(e1_sources.desig_02[i][1:8]))
        e1_b.append(float(e1_sources.desig_02[i][9:]))
        #print(e1_l, e1_b)
        #raise KeyboardInterrupt

    deepglimpse_files=glob.glob('/astro/sargent/bulgeadap/glimpse/deepglimpse/*C*.tbl')

    for i in range(len(deepglimpse_files)):
        print(deepglimpse_files[i])

        deep_sources = read_file(deepglimpse_files[i])

        # convert coordinates to arcseconds
        e1_coords = []
        alle_coords = []
        for i in range(len(e1_l)):
            e1_coords.append([e1_l[i], e1_b[i]])
        for i in range(len(deep_sources.l)):
            alle_coords.append([deep_sources.l.iloc[i], deep_sources.b.iloc[i]])

        # Initialize new all epoch_catalog
        e1_coords = np.transpose(e1_coords)
        deep_coords = np.transpose(deep_coords)

        print(e1_coords)
        raise KeyboardInterrupt
        deep_matches = find_matches(e1_coords, deep_coords)

        #nonoverlap_allepoch = alle_sources.drop(allepoch_matches)

    # write query file
    if os.path.exists(new_cat):
        print('{} and {} already exists, removing and recreating'.format(new_cat))
        os.remove(new_cat)
        os.remove(all_data_cat)

    output = open(new_cat, 'w')
    output.write('\EQUINOX = J2000.0\n')
    output.write('\t\tdesig_orig\t\t\tra\t\tdec\t\t\t\n')

    output_data = open(all_data_cat, 'w')
    output_data.write('\EQUINOX = J2000.0\n')
    output_data.write(
        '\t\tdesig_orig\t\t\tra\t\tdec\t\t\tl\t\tb\t\tmagj\t\tmagk\t\tmagh\t\tmag3\t\tmag4\t\tmag5\t\tmag8\n')

    for i in range(len(nonoverlap_allepoch.l)):
        output.write(
            '{} {}\t{}\t{}\t{}\n'.format(nonoverlap_allepoch.iloc[i].desig_01, nonoverlap_allepoch.iloc[i].desig_02,
                                         nonoverlap_allepoch.iloc[i].ra, nonoverlap_allepoch.iloc[i].dec, 2.0))

        output_data.write(
            '{} {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                nonoverlap_allepoch.iloc[i].desig_01, nonoverlap_allepoch.iloc[i].desig_02,
                nonoverlap_allepoch.iloc[i].ra, nonoverlap_allepoch.iloc[i].dec,
                nonoverlap_allepoch.iloc[i].l, nonoverlap_allepoch.iloc[i].b, nonoverlap_allepoch.iloc[i].magj,
                nonoverlap_allepoch.iloc[i].magh,
                nonoverlap_allepoch.iloc[i].magk, nonoverlap_allepoch.iloc[i].mag3, nonoverlap_allepoch.iloc[i].mag4,
                nonoverlap_allepoch.iloc[i].mag5,
                nonoverlap_allepoch.iloc[i].mag8
            )
        )

    output.close()
    output_data.close()

    print('Saved to {} and {}'.format(new_cat, all_data_cat))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
Description
------------------------------------


Locations
------------------------------------
  The variable save_dir is hardcoded to be /user/ariley/science/ because
  as of this writing (09 October 2017), I am the only person using it.
  If/when this code changes hands, I will try to make that an input.

Arguments
------------------------------------


Examples
------------------------------------

------------------------------------
''' ) )

    parser.add_argument("-g",
                        "--glimpse_data",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name and path of the data from GLIMPSE 1, 2, and 3D to be compared.",
                        default=None,
                        dest = "glimpse_data")

    parser.add_argument("-o",
                        "--output",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name and path of the output file.",
                        default=None,
                        dest="output")



    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    catalog_diff(args.glimpse_data[0], args.output[0])