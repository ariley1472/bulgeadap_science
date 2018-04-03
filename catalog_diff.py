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
#import time

def convert_to_arcsec(long, lat):
    l1 = Angle(long, u.deg)
    b1 = Angle(lat, u.deg)

    return [l1.arcsec, b1.arcsec]

def read_files(epoch1_file, allepoch_file):
    """Read in the files to compare using Pandas DataFrames. These files should be the
    results glimpse3dirolist.py --> *_querylist_alldata.txt"""
    #what is the important information that needs to be read in? Can I just

    allepoch_queryfile = allepoch_file.replace('_alldata.txt', '.txt')
    epoch1_queryfile = epoch1_file.replace('_alldata.txt', '.txt')
    print allepoch_queryfile

    usercols = ['desig_01', 'desig_02', 'l', 'b','ra', 'dec']
    epoch1_sources = read_table(epoch1_file, header = 0, names = usercols,
                               comment = '\\', sep = '\s+',
                               usecols = [0, 1, 2, 3, 4, 5])#12, 13, 14, 15])

    allepoch_sources = read_table(allepoch_file, header = 1, names = usercols,
                                   comment = '\\', sep = '\s+',
                                   usecols = [0, 1, 2, 3, 4, 5])#[3, 4, 12, 13, 14, 15])

    usercols1 = ['desig_01', 'desig_02', 'l', 'b' , 'ra', 'dec', 'magj', 'dmagj', 'magh', 'dmagh', 'magk', 'dmagk', 'mag3', 'dmag3', 'mag4', 'dmag4', 'mag5', 'dmag5', 'mag8', 'dmag8']
    usercols2 = ['desig_01', 'desig_02', 'ra', 'dec', 'major']
    allepoch_diff_df = read_table(allepoch_file, header = 0, names = usercols1,
                              comment = '\\', sep = '\s+',
                              usecols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
                              
    allepoch_df = read_table(allepoch_queryfile, header = 1, names = usercols2, comment = '\\',
                                  sep = '\s+', usecols = [0, 1, 2, 3, 4])
                                  
    print len(allepoch_diff_df)
    print allepoch_diff_df.head(n=10)
    print '-------------------------------------'
    print len(allepoch_df)
    print allepoch_df.head(n = 10)
    print '--------------------------------------'
    
    print len(epoch1_sources)
    print allepoch_diff_df.head(n=10)
    print '-------------------------------------'
    print len(allepoch_sources)
    print allepoch_df.head(n = 10)
    print '--------------------------------------'
    #raise KeyboardInterrupt

    #allepoch_data = readtable()

    #allepoch_querylist = allepoch_file.replace('_querylist_diff.txt', '_querylist.txt')
    
    #allepoch_diff = allepoch_file.replicate('_diff.txt', '_data.txt')
    
    #if not os.path.exists(allepoch_querylist):
        #open(allepoch_querylist, 'a').close()
    #Read in the sources from the IRSA query list
    #lines = []
    #with open(allepoch_querylist, "rw+") as f:
        #for line in f:
            #print line
            #lines.append(line)
            

    #print epoch1_sources
            
    #raise KeyboardInterrupt

    #find the overlap areas on the edges (9 <= b <= 11 & 351 <= b <= 349)
    #print '-------------------------------'
    #print '----FINDING OVERLAP SOURCES----'
    #print '-------------------------------'

    # These criteria are for GLIMPSE2 vs. GLIMPSE1
    #overlap_allepoch = np.where(((allepoch_sources.b < 1.) & (allepoch_sources.b > -1.)) &
    #                               ((((allepoch_sources.l <= 11.) & (allepoch_sources.l >= 9.)))|
    #                            ((allepoch_sources.l <= 351.) & (allepoch_sources.l >= 349.0)) |
    #                             (allepoch_sources.l <= 1.5) | (allepoch_sources.l >= 358.5)))[0]

    #overlap_1epoch = np.where(((epoch1_sources.b < 1.) & (epoch1_sources.b > -1.)) &
    #                          (((epoch1_sources.l <= 11.) & (epoch1_sources.l >= 9.)) |
    #                         ((epoch1_sources.l <= 351.) & (epoch1_sources.l >= 349.0)) |
    #                           (allepoch_sources.l <= 1.5) | (allepoch_sources.l >= 358.5)))[0]

    # For GLIMPSE1 data:
    if 'GLMIC_' in allepoch_file and 'GLMDEEP' not in epoch1_file:
        print 'This is GLIMPSE 1'
        version = 1
        overlap_allepoch = np.where((((allepoch_sources.l >= 9.5) & (allepoch_sources.l <= 10.5)) &
                                    ((allepoch_sources.b >=-1.) & (allepoch_sources.b <= 1.))) |
                                    (((allepoch_sources.b >=-1.) & (allepoch_sources.b <= 1.)) &
                                    ((allepoch_sources.l >= 349.5) & (allepoch_sources.l <= 350.5))))[0]

        #print overlap_allepoch
        overlap_1epoch = np.where((((epoch1_sources.l >= 9.5) & (epoch1_sources.l <= 10.5)) &
                                    ((epoch1_sources.b >=-1.) & (epoch1_sources.b <= 1.))) |
                                    (((epoch1_sources.b >=-1.) & (epoch1_sources.b <= 1.)) &
                                    ((epoch1_sources.l >= 349.5) & (epoch1_sources.l <= 350.5))))[0] #These said "allepoch_sources"
    elif 'v2.0_GLMII_' in allepoch_file:
        print 'This is GLIMPSE 2'
        version = 'GALCEN'
        overlap_allepoch = np.where((((allepoch_sources.l >= 0.) & (allepoch_sources.l <= 1.5)) &
                                    ((allepoch_sources.b >=-1.) & (allepoch_sources.b <= 1.))) |
                                    (((allepoch_sources.b >=-1.) & (allepoch_sources.b <= 1.)) &
                                    ((allepoch_sources.l >= 358.5) & (allepoch_sources.l <= 360.))))[0]
        overlap_1epoch = np.where((((epoch1_sources.l >= 0.) & (epoch1_sources.l <= 1.5)) &
                                    ((epoch1_sources.b >=-1.) & (epoch1_sources.b <= 1.))) |
                                    (((epoch1_sources.b >=-1.) & (epoch1_sources.b <= 1.)) &
                                    ((epoch1_sources.l >= 358.5) & (epoch1_sources.l <= 360.))))[0]

    # Next, the GLIMPSE3D data, but I can't figure out what the swatches on the sky are.. 
    elif 'GLM3D' in allepoch_file:
        print 'This is GLIMPSE 3D'
        version = 'GLIMPSE3D'
        overlap_allepoch = np.where(((allepoch_sources.b >= -4.) & (allepoch_sources.b <= 4.)) &
                                      (((allepoch_sources.l >= 0.) & (allepoch_sources.l <= 33.)) |
                                      ((allepoch_sources.l >= 327.) & (allepoch_sources.l <= 360.))))[0]
        overlap_1epoch = np.where(((epoch1_sources.b >= -4.) & (epoch1_sources.b <= 4.)) &
                                      (((epoch1_sources.l >= 0.) & (epoch1_sources.l <= 33.)) |
                                      ((epoch1_sources.l >= 327.) & (epoch1_sources.l <= 360.))))[0] 
    #print 'allepoch_overlap:', len(overlap_allepoch), min(overlap_allepoch), max(overlap_allepoch)
    print '1epoch_overlap:', overlap_1epoch#len(overlap_1epoch), min(overlap_1epoch), max(overlap_1epoch)

    # Find and count the overlapping sources.
    count = 0
    allepoch_indices = []
    #for i in overlap_1epoch:
        #search_ind = search(epoch1_sources.)

    #for i in overlap_1epoch:
    #    for j in overlap_allepoch:
    #        if #epoch1_sources.desig_02[i] == allepoch_sources.desig_02[j]:#epoch1_sources.desig_02[i] is within 3 arcseconds of allepoch_sources.desig_02[i] or something #epoch1_sources.desig_02[i] == allepoch_sources.desig_02[j]:
    #            print 'matches:', epoch1_sources.desig_02[i], allepoch_sources.desig_02[j]
    #            allepoch_indices.append(j)
    #            count+=1
    #            # drop them from the all-epoch data frame?
    #            # df.drop()? <- does this mean that I need to read in the entire data frame and not
    #            # just the important columns?

    #convert coordinates to arcseconds:
    short_l_epoch1 = epoch1_sources.l[overlap_1epoch]
    short_b_epoch1 = epoch1_sources.b[overlap_1epoch]
    short_ra_epoch1 = epoch1_sources.ra[overlap_1epoch]
    short_dec_epoch1 = epoch1_sources.dec[overlap_1epoch]
    short_desig01_epoch1 = epoch1_sources.desig_01[overlap_1epoch]
    short_desig02_epoch1 = epoch1_sources.desig_02[overlap_1epoch]

    short_l_allepoch = allepoch_sources.l[overlap_allepoch]
    short_b_allepoch = allepoch_sources.b[overlap_allepoch]
    short_ra_allepoch = allepoch_sources.ra[overlap_allepoch]
    short_dec_allepoch = allepoch_sources.dec[overlap_allepoch]
    short_desig01_allepoch = allepoch_sources.desig_01[overlap_allepoch]
    short_desig02_allepoch = allepoch_sources.desig_02[overlap_allepoch]

    print len(overlap_allepoch)
    print len(short_l_allepoch)

    #print len(np.array(lines)[overlap_allepoch])


    #short_lines = np.array(lines)[overlap_allepoch]


    epoch1_coords_short = []
    for i in range(len(short_l_epoch1)):#epoch1_sources.l)):
        epoch1_coords_short.append(convert_to_arcsec(short_l_epoch1.iloc[i], short_b_epoch1.iloc[i]))
    allepoch_coords_short = []
    for i in range(len(short_l_allepoch)):
        allepoch_coords_short.append(convert_to_arcsec(short_l_allepoch.iloc[i], short_b_allepoch.iloc[i]))

    print allepoch_coords_short[0][0], epoch1_coords_short[0][0]
    print allepoch_coords_short[0][1], epoch1_coords_short[0][1]
    print allepoch_coords_short[0], epoch1_coords_short[1]

    print len(allepoch_coords_short)


    #initalize new all epoch catalog file:
    new_cat = allepoch_file.replace('_querylist_alldata.txt', '___querylist_new_test.txt') #formatted to be accepted by IRSA
    bookkeeping = allepoch_file.replace('_querylist_alldata.txt', '___querylist_data_test.txt') #holds ra, dec, l, b
    if os.path.exists(new_cat):
        print 'This file {}  already exists. Removing and Recreating.'.format(new_cat)
        os.remove(new_cat)
        os.remove(bookkeeping)

    #formyinfo = open(bookkeeping, 'w')
    #formyinfo.write('\ EQUINOX = J2000.0\n')
    #formyinfo.write("|     desig           |   glon    |   glat    |   ra  |   dec |\n")

    #f = open(new_cat, 'w')
    #f.write("\ EQUINOX = 'J2000.0'\n")
    #f.write("|         desig           |    ra     |  dec     |  major |\n")
    #f.write("|         char            |   double  |  double  | double |\n")

    # find matches:
    matches_ind = []
    ind = []
    for i in range(len(allepoch_coords_short)):
        print 'i = ', i
        matches = search(allepoch_coords_short[i][0], allepoch_coords_short[i][1], epoch1_coords_short)
        print '------------------'
        print 'matches:', matches

        #matches_ind.append(matches) #index of the matched source in the SINGLE EPOCH CATALOG
        if len(matches) == 0 :
            # This means that there are no matches, so we want to use the information from the all epoch catalog.
            #write data from allepoch_coords[i] out into a new all epoch catalog
            # The following block is to keep information I might find useful for troubleshooting
            #formyinfo.write('{} {}\t{}\t{}\t{}\t{}\n'.format(short_desig01_allepoch.iloc[i], short_desig02_allepoch.iloc[i], #allepoch_sources.desig_01[i], allepoch_sources.desig_02[i],
                                                             #short_l_allepoch.iloc[i], short_b_allepoch.iloc[i], #allepoch_sources.l[i], allepoch_sources.b[i],
                                                             #short_ra_allepoch.iloc[i], short_dec_allepoch.iloc[i]))#allepoch_sources.ra[i], allepoch_sources.dec[i]))
            #Option 2: no matches = leave in all epoch catalog = do nothing.
            print 'No match here!'

            #The following block is to create the file that can be accepted by IRSA query.
            #f.write(short_lines[i])
            #f.write('{} {}  {} {}      {}\n'.format(short_desig01_allepoch.iloc[i], short_desig02_allepoch.iloc[i],
                                                  #short_ra_allepoch.iloc[i], short_dec_allepoch.iloc[i], 2.0))

        else:
            matches_ind.append(matches) #index of the matched source in the SINGLE EPOCH CATALOG
            # This means that there are matches, so we want to use the data from the single epoch catalog.
            #Do not write the data from the allepoch_coords out to the new all epoch catalog.


            #option 2: match = take out of all epoch catalog and leave in single epoch catalog = drop entry from allepoch_sources.
            print 'overlap_allepoch[i]:', overlap_allepoch[i]
            print 'allepoch:', allepoch_sources.l[overlap_allepoch[i]], allepoch_sources.b[overlap_allepoch[i]]
            print 'single epoch:', epoch1_sources.l[matches[0]], epoch1_sources.b[matches[0]]
            print len(allepoch_sources.b), len(overlap_allepoch)
            print 'length before:', len(allepoch_sources)

            ind.append(overlap_allepoch[i])
            #print 'length after:', len(ind)


            print allepoch_coords_short[i][0], allepoch_coords_short[i][1]
            print epoch1_coords_short[matches[0]]
            print 'We have a match!'
            count +=1
    print 'ind:' ,ind
    nonoverlap_allepoch = allepoch_df.drop(ind)
    print len(ind)
    print 'len(allepoch_diff_df):', len(allepoch_diff_df)
    #print 'allepoch_diff_df[17129]:, allepoch_diff_df[17128]
    #raise KeyboardInterrupt
    nonoverlap_allepoch_data =  allepoch_diff_df.drop(ind)
    print 'allepoch_sources shape:', allepoch_sources.shape
    print 'nonoverlap_allepoch shape:', nonoverlap_allepoch.shape
    print 'nonoverlap_allepoch_data shape:', nonoverlap_allepoch_data.shape
    print nonoverlap_allepoch
    
    #print allepoch_sources
    #nonoverlap_allepoch.to_csv(new_cat, header  = "\ EQUINOX = 'J2000.0'\n|         desig           |    ra     |  dec     |  major |\n|         char            |   double  |  double  | double |", index=None, sep='\t')
    #np.savetxt(new_cat, nonoverlap_allepoch.values, fmt=('%s'), header  = "\ EQUINOX = 'J2000.0'\n|         desig           |    ra     |  dec     |  major |\n|         char            |   double  |  double  | double |\n")
    np.savetxt(new_cat, nonoverlap_allepoch.values, fmt=('%s'),
               header="\ EQUINOX = 'J2000.0'\n|         desig           |   ra      |  dec     |  major |\n|         char            |   double  |  double  | double |",
               comments = '', delimiter = '\t')

    np.savetxt(bookkeeping, nonoverlap_allepoch_data.values, fmt = ('%s'), 
               header = "\ EQUINOX = J2000.0\n        desig_orig            gal long    gal lat        ra        dec            mag J    dmag J     mag H    dmag H    mag K    dmag K    mag3    dmag3    mag4    dmag4    mag5    dmag5    mag8    dmag8",
               comments = '', delimiter = '\t')
    print 'saved to {} and {}'.format(new_cat, bookkeeping)
    #formyinfo.close()
    #f.close()

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

        closest_index = [close_indices[np.where(np.array(all_dist) == min(all_dist))[0][0]]]
        print 'closest index:', closest_index, min(all_dist)
        print 'all distances:', all_dist
        #raise KeyboardInterrupt
    else:
        closest_index = close_indices

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

    read_files(args.single_epoch_file[0], args.all_epoch_file[0])