#! /usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from pandas import *
import pandas as pd
import textwrap

def glimpse3dirolist(magnitudes_file, sourcelist, version):
    dir_name, magname = os.path.split(magnitudes_file) #name of magnitude list and the name of the 
    listname = os.path.split(sourcelist)[1] #name of source list
    output = magname.replace('_magnitudes.txt', '_querylist.txt') #normal output # AER 2 Aug
    #output_diff = magname.replace('_magnitudes.txt', '_querylist_diff.txt') # difference output
    output_data = magname.replace('_magnitudes.txt', '_querylist_alldata.txt')
    

    if version == '1':
        plot_title = 'GLIMPSE 1'
    elif version == '2':
        plot_title = 'GLIMPSE II'
    elif version == '3d':
        plot_title = 'GLIMPSE 3D'
    else:
        plot_title = 'DEEP GLIMPSE'

    #Read in the magnitudes:
    #glon, glat, ra, dec = np.loadtxt(magnitudes_file, usecols = (2, 3, 4, 5),
    #                        unpack = True)

    usercols = ['desigorig1', 'desigorig2', 'glon', 'glat', 'ra', 'dec',  'jmag', 'djmag', 'hmag', 'dhmag', 'kmag', 'dkmag',
                'i1mag', 'di1mag', 'i2mag', 'di2mag', 'i3mag', 'di3mag', 'i4mag', 'di4mag']
    mag_sources = read_table(magnitudes_file, header = 0, names = usercols, comment = "\\",
                             sep = '\s+', skiprows = 0,
                             usecols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    print mag_sources.desigorig1[0], mag_sources.desigorig2[0]
    desigorig = []
    for i in range(len(mag_sources.desigorig1)):
        desigorig.append(' '.join([mag_sources.desigorig1[i], mag_sources.desigorig2[i]]))

    #find the red sources
    redind = np.where((mag_sources.i2mag - mag_sources.i4mag) >= 1.0)[0] #find where this condition is true
    numred = len(redind)
    print('numred:', numred)
    

    #plot locations of red sources
    plt.scatter(mag_sources.glon[redind], mag_sources.glat[redind])
    plt.title('Locations of Red Sources')
    plt.savefig(os.path.join(dir_name, '_source_locations.png'))

    # Find possible X-AGBs
    xagbind = np.where(mag_sources.i2mag[redind] <= 7.8)[0] #I think this is finding where the X-AGBs are
    print xagbind
    numxagb = len(xagbind) # Number of X-AGBs
    
    #plot locations of X-AGBs
    plt.scatter(mag_sources.glon[redind[xagbind]],mag_sources.glat[redind[xagbind]])
    plt.title('Locations of Possible X-AGBs')
    plt.savefig(os.path.join(dir_name, '_xagb_locations.png'))
    
    #Read in the sources
    lines = []
    with open(sourcelist, "rw+") as f:
        for line in f:
            if line.startswith('\\') or line.startswith('|'):
                continue
            else:
                lines.append(line)

    print lines[0]
    print lines[1]
    outfile = os.path.join(dir_name, output)
    outfile_data = os.path.join(dir_name, output_data)

    print '{} and {} are the outputs'.format(outfile, outfile_data)

    if os.path.exists(outfile):
        print 'This file {} and {} already exists. Removing and Recreating.'.format(outfile, outfile_data)
        os.remove(outfile)
        os.remove(outfile_data)

    f_data = open(outfile_data, 'w')
    f_data.write('\EQUINOX = J2000.0\n')
    f_data.write('\t\tdesig_orig\t\t\tgal long\tgal lat\t\tra\t\tdec\t\t\tmag J\tdmag J\t mag H\tdmag H\tmag K\tdmag K\tmag3\tdmag3\tmag4\tdmag4\tmag5\tdmag5\tmag8\tdmag8\n')

    #f_diff = open(outfile_diff, 'w')
    #f_diff.write('\ EQUINOX = J2000.0\n')
    #f_diff.write("|     desig           |   glon    |   glat    |   ra  |   dec |\n")

    f = open(outfile, 'w')
    f.write("\ EQUINOX = 'J2000.0'\n")
    f.write("|         desig           |   ra      |  dec     |  major |\n")
    f.write("|         char            |   double  |  double  | double |\n")
    for i in range(len(redind)):
        print '--------------------------------'
        print redind[i]
        print lines[redind[i]]
        print mag_sources.desigorig1[redind[i]], mag_sources.desigorig2[redind[i]]
        f.write(lines[redind[i]]) # AER 2 Aug

        f_data.write('{} {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                mag_sources.desigorig1[redind[i]], mag_sources.desigorig2[redind[i]], mag_sources.glon[redind[i]],
                mag_sources.glat[redind[i]], mag_sources.ra[redind[i]], mag_sources.dec[redind[i]],
                mag_sources.jmag[redind[i]], mag_sources.djmag[redind[i]], mag_sources.hmag[redind[i]],
                mag_sources.dhmag[redind[i]], mag_sources.kmag[redind[i]], mag_sources.dkmag[redind[i]],
                mag_sources.i1mag[redind[i]], mag_sources.di1mag[redind[i]], mag_sources.i2mag[redind[i]],
                mag_sources.di2mag[redind[i]], mag_sources.i3mag[redind[i]], mag_sources.di3mag[redind[i]],
                mag_sources.i4mag[redind[i]], mag_sources.di4mag[redind[i]]))
        print i
    f.close()
    f_data.close()
    print 'saved {} and {}'.format(f, f_data)
    
    print('end')

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

    parser.add_argument("-g", #maybe one of these days, I can have it strip the header for this info... but that's not today.
                        "--vglimpse",
                        action='store',
                        nargs=1,
                        type=str,
                        help="Which glimpse version you're inputting",
                        default=None,
                        dest = "vglimpse")
    parser.add_argument("-m",
                        "--mag_file",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name and path of the magnitude text file to be processed",
                        default=None,
                        dest = "mag_file")
    parser.add_argument("-l",
                        "--list_file",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name and path of the text file containing the list of sources to be processed",
                        default=None,
                        dest = "list_file")

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    print('Files to be processed: {} and {}'.format(args.mag_file[0], args.list_file[0])) #Just checkin'
    print('This file is from GLIMPSE {}'.format(args.vglimpse[0]))

    mag_path, mag_file = os.path.split(args.mag_file[0]) # Breaking up the path and the filename.
    list_path, list_file = os.path.split(args.list_file[0]) #list_path and mag_path should be the same.
    if mag_path != list_path:
        raise IOError('The magnitudes and the list inputs should be from the same path.')


    glimpse3dirolist(args.mag_file[0], args.list_file[0], args.vglimpse[0])
