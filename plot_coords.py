#! /usr/bin/env python

import textwrap
import argparse
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.coordinates import Angle
import astropy.units as u
import numpy as np
from pandas import *


def plot_coords(filename, output):
    print(filename)
    if filename.endswith('.txt'):
    	try:
        	data = ascii.read(filename, data_start=2) #, format='fixed_width', delimiter='\t', comment='#')
        	#print(data[0:10])
        	targs = data['desig']
        	
        except:
        	data = ascii.read(filename, data_start=3)
        	targs = data['col2']
        #targs = targs.s
    elif filename.endswith('.tbl'):
        usercols = ['desig']#orig1', 'desigorig2g']
        targs = read_table(filename, header=3, names=usercols, comment="\\",
                             sep='\s+',
                             usecols=[1])
        #print(data)
        #raise KeyboardInterrupt
        #targs = data['desigorig2']
        #print(targs)
        #raise KeyboardInterrupt
        #data = ascii.read(filename)#, data_start=13)
        #targs = data['designation']

    #print(targs.iloc[0])
    #print('targ', type(targs.iloc[0]))
    #print('type:', type(str(targs.iloc[0]['desig'])))
    #print(targs[0].split(' ')[1][1:])
    #raise KeyboardInterrupt
    long_wrap = []
    lat = []
    for i in range(len(targs)):
    	if filename.endswith('.tbl'):
    		print targs.iloc[i]['desig'][1:]
    		#raise KeyboardInterrupt
        	name = str(targs.iloc[i]['desig'])[1:]#.split(' ')[1][1:]
        else:
        	try:
        		name = targs[i].split(' ')[1][1:]
        	except:
        		name = targs[i].split(' ')[0][1:]
        print('name', name)
        #raise KeyboardInterrupt
        if '-' in name:
            gal_long, gal_lat = name.split('-')
            gal_lat = '-{}'.format(gal_lat)
        elif '+' in name:
            gal_long, gal_lat = name.split('+')
            gal_lat = '+{}'.format(gal_lat)

        gal_long = float(gal_long)
        gal_lat = float(gal_lat)

        long_wrap.append(Angle([gal_long], unit=u.deg).wrap_at('180d').degree)
        lat.append(gal_lat)

    plt.figure(figsize=(15, 3))
    plt.scatter(long_wrap, lat, marker='.', s=0.1)

    plt.ylabel('Galactic Latitude (degrees)')
    plt.xlabel('Galactic Longitude (degrees)')

    plt.savefig(output)#/user/ariley/science/glimpse3d/epoch1_glimpse3d_locations.png')

    print('saved to {}'.format(output))


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
Description
------------------------------------

''' ) )

    parser.add_argument("-f", #maybe one of these days, I can have it strip the header for this info... but that's not today.
                        "--filename",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name of the file containing the list of coordinates to plot",
                        default=None,
                        dest = "filename")
    parser.add_argument("-o", #maybe one of these days, I can have it strip the header for this info... but that's not today.
                        "--output",
                        action='store',
                        nargs=1,
                        type=str,
                        help="output name",
                        default=None,
                        dest = "output")
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    plot_coords(args.filename[0], args.output[0])