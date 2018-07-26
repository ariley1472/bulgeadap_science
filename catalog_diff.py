#! /usr/bin/env python

from astropy import units as u
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
import argparse
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
    usercols = ['desig_01', 'desig_02', 'l', 'b', 'ra', 'dec']
    sources = read_table(filename, header=1, names=usercols, comment='\\', sep='\s+', 
                         usecols=[0, 1, 2, 3, 4, 5])
    
    print(sources)
    return sources
    
def find_matches(epoch1_coords, allepoch_coords):
    #matches_ind = []
    allepoch_matches_ind = []
    ind = []
    count = 0

    epoch1_cat = SkyCoord(epoch1_coords[0], epoch1_coords[1], unit='deg', frame='galactic')
    for i in range(len(allepoch_coords[0])):
        print('i = ', i)
        target = SkyCoord([allepoch_coords[0][i],], [allepoch_coords[1][i],], unit='deg', frame='galactic')
        idxsearcharound, epoch1_match, sep2d, dist3d = epoch1_cat.search_around_sky(target, 2*u.arcsec)

        if len(epoch1_match) == 0: # No matches
            print('No Match Here!')
        elif len(epoch1_match) == 1: # One match
            print('We Have A Match!')
            print(idxsearcharound, epoch1_match, sep2d*u.arcsec)
            print(epoch1_coords[0][epoch1_match], epoch1_coords[1][epoch1_match])
            print(allepoch_coords[0][i], allepoch_coords[1][i])
            print('*******************************************')
            count+=1
            #raise KeyboardInterrupt
            allepoch_matches_ind.append(i)

    print('length:', len(allepoch_matches_ind))

    return allepoch_matches_ind
    

def catalog_diff(epoch1_file, allepoch_file):
    print(type(allepoch_file))
    #raise KeyboardInterrupt

    #read in all the necessary files:
    e1_sources = read_file(epoch1_file)
    alle_sources = read_file(allepoch_file)
    
    #convert coordinates to arcseconds
    e1_coords = []
    alle_coords = []
    for i in range(len(e1_sources.l)):
        e1_coords.append([e1_sources.l.iloc[i], e1_sources.b.iloc[i]])
        #e1_coords.append(convert_to_arcsec(e1_sources.l.iloc[i], e1_sources.b.iloc[i]))
    for i in range(len(alle_sources.l)):
        alle_coords.append([alle_sources.l.iloc[i], alle_sources.b.iloc[i]])
    
    # Initialize new all epoch_catalog
    new_cat = allepoch_file.replace('_querylist_alldata.txt', '_querylist_catdiff.txt')
    if os.path.exists(new_cat): # If it exists, remove so that we can create a separate one.
        os.remove(new_cat)

    e1_coords = np.transpose(e1_coords)
    alle_coords = np.transpose(alle_coords)
    print('shape(e1_coords)', np.shape(e1_coords))

    print(e1_coords)
    #raise KeyboardInterrupt
    
    allepoch_matches = find_matches(e1_coords, alle_coords)

    print('length:', len(allepoch_matches))
    #print(allepoch_matches[0])
    #raise KeyboardInterrupt

    nonoverlap_allepoch = alle_sources.drop(allepoch_matches)

    print(nonoverlap_allepoch.iloc[0].desig_01, nonoverlap_allepoch.iloc[0].desig_02, nonoverlap_allepoch.iloc[0].ra, nonoverlap_allepoch.iloc[0].dec)
    #raise KeyboardInterrupt

    # Create and save these files:
    if os.path.exists(new_cat):
        print('{} already exists, removing and recreating')
        os.remove(new_cat)

    output = open(new_cat, 'w')
    output.write('\EQUINOX = J2000.0\n')
    output.write('\t\tdesig_orig\t\t\tra\t\tdec\t\t\t\n')

    for i in range(len(nonoverlap_allepoch.l)):
        output.write('{} {}\t{}\t{}\t{}\n'.format(nonoverlap_allepoch.iloc[i].desig_01, nonoverlap_allepoch.iloc[i].desig_02,
                                                  nonoverlap_allepoch.iloc[i].ra, nonoverlap_allepoch.iloc[i].dec, 2.0))

    output.close()

    #np.savetxt(new_cat, nonoverlap_allepoch.values, fmt=('%s'),
               #header= "\ EQUINOX = 'J2000.0'\n|         desig           |   ra      |  dec     |  major |\n|         char            |   double  |  double  | double |",
               #comments = '', delimiter = '\t')

    print('Saved to {}'.format(new_cat))


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

    parser.add_argument("-s",
                        "--single_epoch_file",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name and path of the single-epoch catalog to be compared.",
                        default=None,
                        dest = "single_epoch_file")

    parser.add_argument("-a",
                        "--all_epoch_file",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name and path of the all-epoch catalog to be compared.",
                        default=None,
                        dest = "all_epoch_file")



    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    #print args

    catalog_diff(args.single_epoch_file[0], args.all_epoch_file[0])