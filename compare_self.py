#! /usr/bin/env python

#Compare large final_iro_allglimpse.txt to itself to eliminate any redundant entries

from astropy import units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
from astropy.coordinates import ICRS
import numpy as np
import pandas as pd
import sys

def compare_self(cat_file):
    # Read in the catalog using pandas:
    # cat_file = '/user/ariley/science/final_iro_allglimpse.txt'

    usercols = ['desig1', 'desig2', 'ra', 'dec', 'major']
    cat_list = pd.read_table(cat_file, header = 1, names = usercols, comment = "\\",
                             sep = '\s+', usecols = [0, 1, 2, 3, 4])
                             
                             
                             
    print(cat_list.ra)

    cat_data = SkyCoord(ra=cat_list.ra, dec=cat_list.dec, unit='deg', frame = 'icrs')
    matches = np.zeros(len(cat_data.ra), dtype=bool)
    for i in range(len(cat_data.ra)):
    
        c = SkyCoord(cat_data.ra[i], cat_data.dec[i], frame='icrs')
        dist2d = c.separation(cat_data) 
        
        index = np.where(dist2d < 2*u.arcsec)[0]
        dist_2arcsec = dist2d[index].arcsec

        print 'len = {}, {}'.format(len(dist_2arcsec), dist_2arcsec)
        if len(dist_2arcsec) > 2:
            #find the closest one (that is NOT itself):
            min_dist = min(dist_2arcsec[np.where(index != i)[0]])
            closest = index[np.where(dist_2arcsec == min_dist)][0]
            
            print('we have 2 matches, {} matched with {} and {}'.format(i, index[np.where(index != i)[0]][0], index[np.where(index != i)[0]][1]))
            print('{} has a distance away from {}'.format(index[np.where(index != i)[0]][0], dist_2arcsec[np.where(index != i)[0]][0]))
            print('{} has a distance away from {}'.format(index[np.where(index != i)[0]][1], dist_2arcsec[np.where(index != i)[0]][1]))
            print('The closest one is {}, which we will use'.format(closest, min_dist))
            
            matches[closest]
            
        elif len(dist_2arcsec) == 2:
            match_ind = index[np.where(index != i)[0]][0]
            print 'match_ind = ', match_ind
            print('we have a match, {} and {}'.format(i, match_ind))
            print '-----------------------------------'
            if matches[i] == 0 and matches[match_ind] == 0:
                #make them both 1
                matches[i] = 1
                matches[match_ind] = 1
                
            else: # make one of them 0
                   matches[match_ind] = 0
    
    new_cat = cat_list.drop(np.where(matches == True)[0])

    print('length after dropping: {}'.format(len(new_cat.ra)))
    
    output = open('/user/ariley/science/querylist_all_final.txt', 'w')
    output.write("\ EQUINOX = 'J2000.0'\n")
    output.write("|         desig           |   ra      |  dec     |  major |\n")
    output.write("|         char            |   double  |  double  | double |\n")

    for i in range(len(new_cat.ra)):
        row = new_cat.iloc[i]

        output.write(' {} {}  {:10.6f} {:10.6f}      {:2.1f}\n'.format(row.desig1, row.desig2, row.ra, row.dec, row.major))

    # Also read in the full data thing... maybe deal with that later.
    
if __name__ == '__main__':
    cat = sys.argv[1]
    print(cat, type(cat))
    compare_self(cat)