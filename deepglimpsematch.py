#! /usr/bin/env python

from astropy import units as u
from astropy.coordinates import Angle
import argparse
import numpy as np
from math import sqrt
#import matplotlib.pyplot as plt
import os
from pandas import *
import pandas as pd
import scipy.spatial as spatial
import textwrap

def convert_to_arcsec(long, lat):
    l1 = Angle(long, u.deg)
    b1 = Angle(lat, u.deg)

    return [l1.arcsec, b1.arcsec]


def match(other_data, deepglimpse):
    # Read in the GLIMPSE 1 data:
    
    print 'other_data:', other_data
    print 'deepglimpse:', deepglimpse
    deepusercols = ['desigorig1', 'desigorig2', 'glon', 'glat', 'ra', 'dec',
                    'magj', 'dmagj', 'magh', 'dmagh', 'magk', 'dmagk', 'mag3', 'dmag3',
                    'mag4', 'dmag4', 'mag5', 'dmag5', 'mag8', 'dmag8']
                    
    gdeepsources = read_table(deepglimpse, header=3, names=deepusercols, comment="\\", sep='\s+', skiprows=0,
                           usecols=[0, 1, 4, 5, 8, 9, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26])
    
    usercols = ['desigorig1', 'desigorig2', 'glon', 'glat', 'ra', 'dec', 'jmag', 'djmag', 'hmag', 'dhmag', 'kmag',
                'dkmag', 'i1mag', 'di1mag', 'i2mag', 'di2mag', 'i3mag', 'di3mag', 'i4mag', 'di4mag']
    g123dsources = read_table(other_data, header=0, names=usercols, comment="\\", sep='\s+', skiprows=0,
                           usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])

    # Convert to arcseconds
    g123dcoords = []
    for i in range(len(g123dsources.desigorig2)):
        g123dcoords.append(convert_to_arcsec(g123dsources.glon[i], g123dsources.glat[i]))
    gdeepcoords = []
    for i in range(len(gdeepsources.desigorig2)):
    	print gdeepsources.glon[i], gdeepsources.glat[i]
        gdeepcoords.append(convert_to_arcsec(gdeepsources.glon[i], gdeepsources.glat[i]))

    # Initialize the file
    output = deepglimpse.replace('.tbl', '_matches.txt')
    print output
    print deepglimpse
    #raise KeyboardInterrupt
    if not os.path.exists(output):
        out = open(output, 'a')
    else:
        out = open(output, 'w')

    # Write header
    out.write(
        '\t\tdesig_orig\t\t\tgal long\tgal lat\t\tra\t\tdec\t\t\tmag J\tdmag J\t mag H\tdmag H\tmag K\tdmag K\tmag3\tdmag3\tmag4\tdmag4\tmag5\tdmag5\tmag8\tdmag8\n')

    # find matches:
    matches_ind = []
    ind = []
    count = 0
    for i in range(len(g123dcoords)):
        print 'i = ', i
        #print 'Coordinates:', g123dcoords[i][0], g123dcoords[i][1]
        #print 'coords:', gdeepcoords
        #raise KeyboardInterrupt
        matches = search(g123dcoords[i][0], g123dcoords[i][1], gdeepcoords) #find matches

        if len(matches) == 0 :
            print 'No match here!'

        else:
            matches_ind.append(matches)
            # This means that there are matches, so we want to use the data from the single epoch catalog.
            #Do not write the data from the allepoch_coords out to the new all epoch catalog.


            ind.append(i)

            print gdeepcoords[matches[0]][0], gdeepcoords[matches[0]][1]
            print g123dcoords[i]
            print 'length of matches:', len(matches)
            print 'matches:', matches
            print 'type of matches:', type(matches)
            print 'len(g123dcoords):', len(g123dcoords)
            print 'We have a match!'
            count +=1
    out.close()
    print 'saved to {}'.format(new_cat)

    print "We found a total of {} matching sources!".format(count)

def search(xcoord, ycoord, epoch1_coords):#ras, decs):
    '''Finds all the sources from the all epoch catalog that are within 3 arcseconds from (xcoord, ycoord), which are
    coordinates from the single epoch catalog.'''
    print xcoord, ycoord, len(epoch1_coords)

    point_tree = spatial.cKDTree(epoch1_coords) #may need to convert these into arcseconds?
    close_indices = point_tree.query_ball_point([xcoord, ycoord], 2) # 3 arcseconds?
    close_points = point_tree.data[point_tree.query_ball_point([xcoord, ycoord], 2)]

    #Get the *closest* point
    if len(close_points) > 1:
        print 'len(close_points):', len(close_points)
        all_dist = []
        for c in close_points:
            print 'epoch1_coords:', c[0], c[1]
            print 'allepoch_coords:', xcoord, ycoord
            distance = sqrt(abs(c[0] - xcoord)**2 + abs(c[1] - ycoord)**2)
            all_dist.append(distance)

        closest_index = close_indices[np.where(np.array(all_dist) == min(all_dist))[0][0]]
        print 'closest index:', closest_index, min(all_dist)
        print 'all distances:', all_dist
        #raise KeyboardInterrupt
    else:
        closest_index = close_indices

	print 'closest_index:', type(closest_index), closest_index
    return closest_index

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

    parser.add_argument("-d",
                        "--deepglimpse_file",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name and path of the deepglimpse table to be compared.",
                        default=None,
                        dest = "deepglimpse_file")



    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    match(args.glimpse_data[0], args.deepglimpse_file[0])