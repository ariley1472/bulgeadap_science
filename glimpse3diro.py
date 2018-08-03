#! /usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from pandas import *
import pandas as pd
import textwrap
import time

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
        return i + 1

def glimpse3diro(filename_, version, output_dir, half, filetype):
    catname = os.path.split(filename_)[1] #name of catalog (or archive)
    print 'half', half

    if half is not None:
        tot_rows = file_len(filename_)
        if half[0] == 'last':
            skip_rows = tot_rows/2
            n_rows = None
        if half[0] == 'first':
            skip_rows = 0
            n_rows = tot_rows/2
            print 'n_rows = {}'.format(n_rows)

    else:
        print('else')
        skip_rows = 0
        n_rows = None
        print 'n_rows = {}'.format(n_rows)


    if version == '1':
        plot_title = 'GLIMPSE 1'
    elif version == '2':
        plot_title = 'GLIMPSE II'
    elif version == '3d':
        plot_title = 'GLIMPSE 3D'
    else:
        plot_title = 'DEEP GLIMPSE'

    #Read in the lines, using PANDAS! HUZZAH!:
    if filetype == "archive":
        print('File type is archive')
        usercols = ['desigorig1', 'desigorig2', 'glonorig', 'glatorig', 'raorig', 'decorig', 'csforig',
                    'magjorig', 'dmagjorig', 'maghorig', 'dmaghorig', 'magkorig', 'dmagkorig',
                    'mag3orig', 'dmag3orig', 'mag4orig', 'dmag4orig', 'mag5orig', 'dmag5orig',
                    'mag8orig', 'dmag8orig', 'f45orig', 'df45orig', 'f80orig', 'df80orig']
        sources = read_table(filename_, header = 3, names = usercols, comment = "\\",
                             sep = '\s+', nrows = n_rows, skiprows = skip_rows,
                             usecols = [0, 1, 4, 5, 8, 9, 12, 13, 14, 15, 16, 17, 18, 19,
                                        20, 21, 22, 23, 24, 25, 26, 35, 36, 39, 40])
        #iro3dind = np.where((((((((((sources.f45orig >= 0.5) & (sources.df45orig > -900.)) &
                            #(sources.f80orig >= 10.)) & (sources.df80orig > -900.)) &
                            #(sources.f45orig <= 450)) & (sources.f80orig <= 1590)) &
                            #(sources.csforig == 0)) & ((sources.df45orig/sources.f45orig) <= 0.15)) & # changed from <1.5 to ==0 (csforig)
                            #((sources.df80orig/sources.f80orig) <= 0.15)))
        iro3dind = np.where((sources.f45orig >= 0.5) & (sources.df45orig > -900.) &
                            (sources.f80orig >= 10.) & (sources.df80orig > -900.) &
                            (sources.f45orig <= 450.) & (sources.f80orig <= 1590.) &
                            (sources.csforig == 0) & ((sources.df45orig/sources.f45orig) <= 0.15) &
                            ((sources.df80orig/sources.f80orig) <= 0.15))[0]                    



    #print sources.desigorig1[0]

    elif filetype == "catalog":
        print('filetype is catalog')
        print('skiprows = ', skip_rows)
        usercols = ['desigorig1', 'desigorig2', 'glonorig', 'glatorig', 'raorig', 'decorig', 'csforig',
                    'magjorig', 'dmagjorig', 'maghorig', 'dmaghorig', 'magkorig', 'dmagkorig',
                    'mag3orig', 'dmag3orig', 'mag4orig', 'dmag4orig', 'mag5orig', 'dmag5orig',
                    'mag8orig', 'dmag8orig', 'f45orig', 'df45orig', 'f80orig', 'df80orig', 'rms45orig',
                    'rms80orig', 'numdet45', 'numdet80']
        sources = read_table(filename_, header = 3, names = usercols, comment = "\\",
                             sep = '\s+', nrows = n_rows, skiprows = skip_rows,
                             usecols = [0, 1, 4, 5, 8, 9, 12, 13, 14, 15, 16, 17, 18, 19,
                                        20, 21, 22, 23, 24, 25, 26, 35, 36, 39, 40, 42, 44, 61, 63])
        print sources.desigorig1[0]
        #iro3dind = np.where((((((((((((((sources.f45orig >= 0.5) & (sources.df45orig > -900.)) &
        #                    (sources.f80orig >= 10.)) & (sources.df80orig > -900.)) &
        #                    (sources.f45orig <= 450.)) & (sources.f80orig <= 1590.)) &
        #                    (sources.csforig == 0)) & ((sources.df45orig/sources.f45orig) <= 0.15)) & # changed from <1.5 to ==0 (csforig)
        #                    ((sources.df80orig/sources.f80orig) <= 0.15)) & (sources.numdet45 >= 2.0)) &
        #                    (sources.numdet80 >= 2.0)) &
        #                    ((sources.rms45orig/sources.f45orig) <= 0.15)) &
        #                   ((sources.rms80orig/sources.f80orig) <= 0.15)))[0]
        iro3dind = np.where((sources.f45orig >= 0.5) & (sources.df45orig > -900.) &
                            (sources.f80orig >= 10.) & (sources.df80orig > -900.) &
                            (sources.f45orig <= 450.) & (sources.f80orig <= 1590.) &
                            (sources.csforig == 0) & ((sources.df45orig/sources.f45orig) <= 0.15) &
                            ((sources.df80orig/sources.f80orig) <= 0.15) & (sources.numdet45 >= 2.0) &
                            (sources.numdet80 >= 2.0) & ((sources.rms45orig/sources.f45orig) <= 0.15) &
                            ((sources.rms80orig/sources.f80orig) <= 0.15))[0]

    else:
        raise IOError("Please specify the file type and try again.")


    # Plot the locations of the objects and save it or else I think the rest of it will hang:
    plt.scatter(sources.glonorig, sources.glatorig)
    plt.title('{} {}'.format(plot_title, filetype))
    plt.savefig(os.path.join(output_dir, catname.replace('.tbl', '_test_plot1.png')))

    archinput_file = os.path.join(output_dir, catname.replace('.tbl', '_test_inputlist.txt'))
    if os.path.exists(archinput_file):
        print 'This file already exists. Removing and Recreating.'
        os.remove(archinput_file)
    f = open(archinput_file, 'w')
    f.write("\ EQUINOX = 'J2000.0'\n")
    f.write("|         desig           |   ra      |  dec     |  major |\n")
    f.write("|         char            |   double  |  double  | double |\n")
    for i in iro3dind:
        f.write(' {} {}  {:10.6f} {:10.6f}      {:2.1f}\n'.format(sources.desigorig1[i], sources.desigorig2[i], sources.raorig[i], sources.decorig[i], 2.0))
        print i

    f.close()

    #write the archive magnitudes list
    archmags_file = os.path.join(output_dir, catname.replace('.tbl', '_test_magnitudes.txt'))
    if os.path.exists(archmags_file):
        print 'This file already exists. Removing and Recreating.'
        os.remove(archmags_file)
    f = open(archmags_file, 'w')
    f.write('\t\tdesig_orig\t\t\tgal long\tgal lat\t\tra\t\tdec\t\t\tmag J\tdmag J\t mag H\tdmag H\tmag K\tdmag K\tmag3\tdmag3\tmag4\tdmag4\tmag5\tdmag5\tmag8\tdmag8\n')
    n=0
    for i in iro3dind: #range(len(iro3dind)):
        print n
        f.write('{} {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                sources.desigorig1[i], sources.desigorig2[i], sources.glonorig[i], sources.glatorig[i], sources.raorig[i],
                sources.decorig[i], sources.magjorig[i], sources.dmagjorig[i], sources.maghorig[i],
                sources.dmaghorig[i], sources.magkorig[i], sources.dmagkorig[i], sources.mag3orig[i],
                sources.dmag3orig[i], sources.mag4orig[i], sources.dmag4orig[i], sources.mag5orig[i],
                sources.dmag5orig[i], sources.mag8orig[i], sources.dmag8orig[i]))
        n+=1
        #print i

    f.close()

    return 'end'

def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=textwrap.dedent('''\
Description
------------------------------------
  Automated script to go through an input GLIMPSE catalog or archive.
  It reads in any necessary columns in the input table and finds
  which sources in the table satisfy the Robitaille, et al 2008
  criteria. Then, it outputs the magnitudes and a list of the sources
  in magitude and list text files, respectively, to be used as input
  for the next script, glimpse3dirolist.p(ro or y). The python version
  is in development as of 09 October 2017.

Locations
------------------------------------
  The variable save_dir is hardcoded to be /user/ariley/science/ because
  as of this writing (09 October 2017), I am the only person using it.
  If/when this code changes hands, I will try to make that an input.

Arguments
------------------------------------
-f or --filename: the archive or catalog to be processed
-g or --vglimpse: The glimpse version (options are 1, 2, 3d, & deep)

Examples
------------------------------------
  The following command tells glimpse3diro.py to run on the file
  /astro/sargent/bulgeadap/glimpse/deepglimpse/GLMDEEPC_l024-026.tbl, which is
  a DEEP GLIMPSE archive.

  $ python glimpse3diro.py -f /astro/sargent/bulgeadap/glimpse/deepglimpse/GLMDEEPC_l024-026.tbl -v deep

  The following command tells glimpse3diro.py to run on the file
  /astro/sargent/bulgeadap/glimpse/glimpse3d/tiles/GLM3DA_1335-02.tbl, which is
  a GLIMPSE 3D archive.

  $ python glimpse3diro.py -f /astro/sargent/bulgeadap/glimpse/glimpse3d/tiles/GLM3DA_1335-02.tbl -g 3d
------------------------------------
''' ) )

    parser.add_argument("-v", #maybe one of these days, I can have it strip the header for this info... but that's not today.
                        "--vglimpse",
                        action='store',
                        nargs=1,
                        type=str,
                        help="Which glimpse version you're inputting",
                        default=None,
                        dest = "vglimpse")
    parser.add_argument("-f",
                        "--filename",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name and path of the archive or catalog to be processed",
                        default=None,
                        dest = "filename")
    parser.add_argument("-d",
                        "--do_half",
                        action = 'store',
                        nargs = 1,
                        type = str,
                        help = "If the archive or catalog is too large, consider running the first half and then the second half",
                        default = None,
                        dest = "do_half")
    parser.add_argument("-t",
                        "--filetype",
                        nargs=1,
                        type=str,
                        help="Is this an archive or catalog?",
                        default=None,
                        dest="filetype")

    return parser.parse_args()

if __name__ == '__main__':
    t = time.time()
    args = parse_args()
    print(args)
    print(args.filetype[0])
    #raise KeyboardInterrupt

    print('File to be processed: {}'.format(args.filename[0])) #Just checkin'
    print('This file is from GLIMPSE {}'.format(args.vglimpse[0]))

    path, filename_short = os.path.split(args.filename[0]) # Just the name of the file.. no path.
    dir_name = filename_short.replace('.tbl', '')
    save_dir = os.path.join('/user/ariley/science/glimpse{}_test'.format(args.vglimpse[0]), dir_name)

    if not os.path.exists(save_dir): #set up the directory structure
        print('making directory {}'.format(save_dir))
        os.mkdir(save_dir)


    print('The outputs will be saved in {}'.format(save_dir))
    glimpse3diro(args.filename[0], args.vglimpse[0], save_dir, args.do_half, args.filetype[0])
    print('This took {} seconds to run'.format(time.time() - t))
