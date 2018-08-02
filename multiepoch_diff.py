#! /usr/bin/env python

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
import argparse
import numpy as np
from math import sqrt
import os
import pandas as pd
from pandas import *
import scipy.spatial as spatial
import textwrap


def multi_epoch_diff(epoch1file, epoch2file, output_file):
    #the input should be querylist.txt or querylist_new.txt, depending upon the specific file. 

    #construct the names used for the output:
    #print epoch1file
    file_split = '/'.join(epoch1file.split('/')[:-2])
    qlist = '{}/{}_full_querylist_logicalor.txt'.format(file_split, output_file)
    dlist = '{}/{}_multiepoch_data.txt'.format(file_split, output_file)
    #nlist = file_split + '/multiepoch_diff_nomatch.txt'
    
    smdata = '{}/{}_full_datalist.txt'.format(file_split, output_file)
    
    #Get any other necessary files
    #datafile_epoch1 = get_other_files(epoch1file)
    #datafile_epoch2 = get_other_files(epoch2file)
        
    #read in the querylist files
    #epoch1_querylist = read_file(epoch1file)
    #e1_query_list = epoch1_querylist.to_string(index = False).split('\n')[1:]
    
    #epoch2_querylist = read_file(epoch2file)
    #e2_query_list = epoch2_querylist.to_string(index = False).split('\n')[1:]
    
    #read in the data files
    #epoch1_alldata = read_file(datafile_epoch1)
    #e1_alldata_list = epoch1_alldata.to_string(index = False).split('\n')[1:]
    
    #epoch2_alldata = read_file(datafile_epoch2)
    #e2_alldata_list = epoch2_alldata.to_string(index = False).split('\n')[1:]
    
    #match them!
    match(epoch1file, epoch2file, qlist, dlist, smdata)#epoch1_alldata, epoch2_alldata, qlist, dlist) #querylist, epoch2_querylist) #, output_files = [qlist, dlist])
    
    return
    
def get_other_files(input_file):

    try:
        data_file = input_file.replace('.txt', '_alldata.txt')
    except:
        data_file = input_file.replace('.txt', '_data.txt')
        
    return data_file
    
def read_file(file):
    
    if file.endswith('_querylist.txt') or file.endswith('_querylist_new.txt'):
        usercols = ['desig_01', 'desig_02', 'ra', 'dec', 'major'] #this may not be right.
        columns = [0, 1, 2, 3, 4]
        
    else:
        usercols = ['desig_01', 'desig_02', 'l', 'b' , 'ra', 'dec', 'magj', 'dmagj', 'magh', 'dmagh', 'magk', 'dmagk', 'mag3', 'dmag3', 'mag4', 'dmag4', 'mag5', 'dmag5', 'mag8', 'dmag8']
        columns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        
    sources = pd.read_table(file, header = 0, names = usercols,
                            comment = '\\', sep = '\s+',
                            usecols = columns)           
    return sources

def match(file1, file2, querylist_name, datalist_name, smalldata_name):#, output_files = ['output.txt']):

    #Get any other necessary files
    datafile_epoch1 = get_other_files(file1)
    datafile_epoch2 = get_other_files(file2)

    #epoch1_querylist = read_file(file1)
    #e1_query_list = epoch1_querylist.to_string(index = False).split('\n')[2:]
    #change to readlines!! 
    e1_query = open(file1)
    e1_query_list = e1_query.readlines()[3:]
    e1_query.close()
    
    #epoch2_querylist = read_file(file2)
    #e2_query_list = epoch2_querylist.to_string(index = False).split('\n')[2:]
    e2_query = open(file2)
    e2_query_list = e2_query.readlines()[3:]
    e2_query.close()
    
    
    #read in the data files
    epoch1_alldata = read_file(datafile_epoch1)
    e1_alldata_list = epoch1_alldata.to_string(index = False).split('\n')[1:]
    ##e1_data = open(datafile_epoch1)
    ##e1_alldata_list = e1_data.readlines()[2:]
    ##e1_data.close()
    
    epoch2_alldata = read_file(datafile_epoch2)
    e2_alldata_list = epoch2_alldata.to_string(index = False).split('\n')[1:]
    ##e2_data = open(datafile_epoch1)
    ##e2_alldata_list = e2_data.readlines()[2:]
    ##e2_data.close()
    
    #print e1_alldata_list[0]
    #print e2_alldata_list[0]
    #print e1_query_list[0]
    #print e2_query_list[0]
    #raise KeyboardInterrupt
    
    cat = SkyCoord(epoch1_alldata.l, epoch1_alldata.b, unit = 'deg', frame = 'galactic') #The catalog we're comparing to is the epoch1 catalog
    #initialize mask:
    mask = np.zeros(len(cat), dtype = bool)
    
    #initialize files:
    if os.path.exists(querylist_name):
        print('Removing {}'.format(querylist_name))
        os.remove(querylist_name)
    if os.path.exists(datalist_name):
        print('Removing {}'.format(datalist_name))
        os.remove(datalist_name)
    #if os.path.exists(nomatch_name):
    #    print('Removing {}'.format(nomatch_name))
    #    os.remove(nomatch_name)
    
    print('Recreating {} and {}'.format(datalist_name, querylist_name))
    query_output = open(querylist_name, 'w')
    query_output.write("\ EQUINOX = 'J2000.0'\n")
    query_output.write("|         desig           |   ra      |  dec     |  major |\n")
    query_output.write("|         char            |   double  |  double  | double |\n")
        
    data_output = open(datalist_name, 'w')
    data_output.write("----------------------------------------------------------------------------EPOCH 1----------------------------------------------------------------------------------------------------      ----------------------------------------------------------------------------EPOCH 2---------------------------------------------------------------------\n")
    #data_output.write("\t\tdesig_orig\t\t\tgal long\tgal lat\t\t\tra\t\tdec\t\tmag J\tdmag J\t mag H\tdmag H\tmag K\tdmag K\tmag3\tdmag3\tmag4\tdmag4\tmag5\tdmag5\tmag8\tdmag8\t||| \tdesig_orig\t\t\t\tgal long\tgal lat\t\tra\t\tdec\t\t\tmag J\tdmag J\t mag H\tdmag H\tmag K\tdmag K\tmag3\tdmag3\tmag4\tdmag4\tmag5\tdmag5\tmag8\tdmag8\n")
        
        
    more_output = open(smalldata_name, 'w')
    more_output.write('desig1	desig2		long    lat     ra		dec		magj	magh	magk	mag3	mag4	mag5	mag8	blah	blah	blah	blah 	BL\n')
    #more_output.write('\t\tname\t\tglon\tglat\tra\t\tdec\t\tmagj\tmagh\tmagk\tmag36\tmag45\tmag80\tmage\tmag24\tfname\n')
    
    #nomatch_output = open(nomatch_name, 'w')
    
    
    n = 0
    for i in range(len(epoch2_alldata.b)):
        print i, epoch2_alldata.l[i], epoch2_alldata.b[i]
        target = SkyCoord([epoch2_alldata.l[i],], [epoch2_alldata.b[i],], unit = 'deg', frame = 'galactic')
        idxsearcharound, epoch1_match, sep2d, dist3d = cat.search_around_sky(target, 2*u.arcsec) #searches for matches within 2 arcseconds away
        print(epoch1_match)
        #print(epoch1_match[0])
        #print epoch1_alldata.desig_01[epoch1_match[0]]
        #raise KeyboardInterrupt
        #print epoch2_alldata.desig_01[1]
        #raise KeyboardInterrupt
        if len(epoch1_match) == 0:
            print('No match found for {}'.format(target))
            
            #append to a list or a file
            #print e2_query_list[i]
            #raise KeyboardInterrupt
            query_output.write('{}'.format(e2_query_list[i]))#e2_)#' {}\n'.format(e2_query_list[i]))#, '\n')
            
            #print epoch2_alldata.desig_01[i]
            
            #also append to some sort of housekeeping file (keeps track of ALL data and inputs -9999 for any data without a match). 
            data_output.write('{}  {} \t\tN/A\t\tN/A\t\t\tN/A\t\t\tN/A\t\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t\t{}\n'.format(epoch2_alldata.desig_01[i], epoch2_alldata.desig_02[i], e2_alldata_list[i]))#data2.desig_01[i], data2.desig_02[i], e2_alldata_list[i]))#tmag J\tdmag J\t mag H\tdmag H\tmag K\tdmag K\tmag3\tdmag3\tmag4\tdmag4\tmag5\tdmag5\tmag8\tdmag8\n')
            more_output.write('{} {}\t{:9f}\t{:10f}\t{:10f}\t{:10f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{}\t{}\t{}\t{}\t{}\n'.format(epoch2_alldata.desig_01[i], epoch2_alldata.desig_02[i],
                                                                                                  epoch2_alldata.l[i], epoch2_alldata.b[i],
                                                                                                  epoch2_alldata.ra[i], epoch2_alldata.dec[i], 
                                                                                                  epoch2_alldata.magj[i], epoch2_alldata.magh[i], epoch2_alldata.magk[i],
                                                                                                  epoch2_alldata.mag3[i], epoch2_alldata.mag4[i], epoch2_alldata.mag5[i],
                                                                                                  epoch2_alldata.mag8[i], 
                                                                                                  '       ', '       ', '       ', '       ', 'XX')) 
            #nomatch_output.write('{}'.format(e2_query_list[i]))
        elif len(epoch1_match) == 1:
            print('1 match found for {}'.format(target))
            #print epoch1_match[0]
            
            #append to a list or a file
            query_output.write('{}'.format(e2_query_list[i]))#, '\n')
            
            #also append to some sort of housekeeping file (keeps track of ALL data and inputs -9999 for any data without a match).
            data_output.write('{}\t\t{}\n'.format(e1_alldata_list[epoch1_match], e2_alldata_list[i]))#{} {}\t{}\t{}\t{}\t{}\t-99999\t-99999\t -99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t\t \t\t {}\n'.format(epoch1_alldata.desig_01[epoch1_match[0]], epoch1_alldata.desig_02[epoch1_match[0]], epoch1_alldata.l[epoch1_match[0]], epoch1_alldata.b[epoch1_match[0]], epoch1_alldata.ra[epoch1_match[0]], epoch1_alldata.dec[epoch1_match[0]], e2_alldata_list[i]))#tmag J\tdmag J\t mag H\tdmag H\tmag K\tdmag K\tmag3\tdmag3\tmag4\tdmag4\tmag5\tdmag5\tmag8\tdmag8\n')
            
            #print epoch1_alldata.desig_01[epoch1_match[0]]
            #print epoch1_alldata.desig_02[epoch1_match[0]], epoch1_alldata.b[epoch1_match[0]], epoch1_alldata.l[epoch1_match[0]]
            #raise KeyboardInterrupt  
            
            more_output.write('{} {}\t{:9f}\t{:10f}\t{:10f}\t{:10f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{}\t{}\t{}\t{}\t{}\n'.format(epoch1_alldata.desig_01[epoch1_match[0]], epoch1_alldata.desig_02[epoch1_match[0]], 
                                                                                                  epoch1_alldata.l[epoch1_match[0]], epoch1_alldata.b[epoch1_match[0]],
                                                                                                  epoch1_alldata.ra[epoch1_match[0]], epoch1_alldata.dec[epoch1_match[0]], 
                                                                                                  epoch1_alldata.magj[epoch1_match[0]], epoch1_alldata.magh[epoch1_match[0]], epoch1_alldata.magk[epoch1_match[0]],
                                                                                                  epoch1_alldata.mag3[epoch1_match[0]], epoch1_alldata.mag4[epoch1_match[0]], epoch1_alldata.mag5[epoch1_match[0]],
                                                                                                  epoch1_alldata.mag8[epoch1_match[0]],
                                                                                                  '       ', '       ', '       ', '       ', 'XX')) 
        
            #start constructing a mask for the catalog
            ## True in the index with a match
            ## False in the index without a match
            mask[epoch1_match[0]] = True
            n+=1
            
        elif len(epoch1_match) > 1:
            print('{} matches found for {}'.format(len(epoch1_match), target))
            
            #find the closest match
            closest_ind = np.where(sep2d == min(sep2d))[0][0]
            #epoch1_ind = epoch2_match[closest_ind] #13 July 2018: commented this out.
            epoch1_ind = epoch1_match[closest_ind] #13 July 2018: added this instead.
            #epoch1_match_closest = cat[epoch1_ind]
            
            #append to a list or a file
            query_output.write('{}'.format(e2_query_list[epoch1_ind]))#, '\n')
            print epoch1_ind
            raise KeyboardInterrupt
            
            #also append to some sort of housekeeping file (keeps track of ALL data and inputs -9999 for any data without a match). 
            data_output.write('{}\t\t{}\n'.format(e1_alldata_list[epoch1_ind], e2_alldata_list[i]))#{} {}\t{}\t{}\t{}\t{}\t-99999\t-99999\t -99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t-99999\t\t \t\t {}\n'.format(epoch1_data.desig_01[epoch1_ind[0]], epoch1_data.desig_02[epoch1_ind[0]], epoch1_data.l[epoch1_ind[0]], epoch1_data.b[epoch1_ind[0]], epoch1_data.ra[epoch1_ind[0]], epoch1_data.dec[epoch1_ind[0]], e2_alldata_list[i]))#tmag J\tdmag J\t mag H\tdmag H\tmag K\tdmag K\tmag3\tdmag3\tmag4\tdmag4\tmag5\tdmag5\tmag8\tdmag8\n'
            
            
            more_output.write('{} {}\t{:9f}\t{:10f}\t{:10f}\t{:10f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{}\t{}\t{}\t{}\t{}\n'.format(epoch1_alldata.desig_01[epoch1_ind], epoch1_alldata.desig_02[epoch1_ind],
                                                                                                  epoch1_alldata.l[epoch1_ind], epoch1_alldata.b[epoch1_ind],
                                                                                                  epoch1_alldata.ra[epoch1_ind], epoch1_alldata.dec[epoch1_ind], 
                                                                                                  epoch1_alldata.magj[epoch1_ind], epoch1_alldata.magh[epoch1_ind], epoch1_alldata.magk[epoch1_ind],
                                                                                                  epoch1_alldata.mag3[epoch1_ind], epoch1_alldata.mag4[epoch1_ind], epoch1_alldata.mag5[epoch1_ind],
                                                                                                  epoch1_alldata.mag8[epoch1_ind],
                                                                                                  '       ', '       ', '       ', '       ','XX')) 
            
            #start constructing a mask for the catalog
            ## True in the index with a match
            ## False in the index without a match
            #mask[epoch2_ind] = True #13 July 2018 commented this out
            mask[epoch1_ind] = True
            n+=1
            
    #To account for the epoch1 data that was NOT a match, but should still be included, 
    #use the mask I created. Anywhere it's FALSE, append that data to a list or a file.
    
    print mask
    
    epoch1_notmatch = np.where(mask == False)[0]
    print len(epoch1_notmatch)
    for i in epoch1_notmatch:
        print 'i in mask:', i, e1_query_list[i]
        #append to a list or a file:
        query_output.write('{}'.format(e1_query_list[i]))#, '\n')
        #something with data1.desig_01[i], data1.desig_02[i], data1.l[i]
        data_output.write('{}\t\t N/A \t\t\t\t\t\tN/A\t\t\tN/A\t\t\tN/A\t\t\tN/A\t\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t99.999\t\t\n'.format(e1_alldata_list[i]))#, data2.desig_01[epoch2], data2.desig_02[epoch2], data2.l[epoch2], data2.b[epoch2], data2.ra[epoch2], data2.dec[epoch2]))
        
        #more_output.write(
        more_output.write('{} {}\t{:9f}\t{:10f}\t{:10f}\t{:10f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:3.3f}\t{:7}\t{:7}\t{:7}\t{:7}\t{:5}\n'.format(epoch1_alldata.desig_01[i], epoch1_alldata.desig_02[i], 
                                                                                                  epoch1_alldata.l[i], epoch1_alldata.b[i],
                                                                                                  epoch1_alldata.ra[i], epoch1_alldata.dec[i], 
                                                                                                  epoch1_alldata.magj[i], epoch1_alldata.magh[i], epoch1_alldata.magk[i],
                                                                                                  epoch1_alldata.mag3[i], epoch1_alldata.mag4[i], epoch1_alldata.mag5[i],
                                                                                                  epoch1_alldata.mag8[i],
                                                                                                  '       ', '       ', '       ', '       ', 'XX')) 
    
    #matches_mask = np.where(mask == False)[0]
    #for i in matches_mask:
        #nomatch_output.write('{}'.format(e1_query_list[i]))
        
    #save_file(filename)
    data_output.close()
    query_output.close()
    #nomatch_output.close()
    
    print('{} matches found!'.format(n))
    return
    

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

    parser.add_argument("-e",
                        "--epoch1",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name and path of the single-epoch catalog to be compared.",
                        default=None,
                        dest = "single_epoch_file")

    parser.add_argument("-t",
                        "--epoch2",
                        action='store',
                        nargs=1,
                        type=str,
                        help="The name and path of the all-epoch catalog to be compared.",
                        default=None,
                        dest = "all_epoch_file")

    parser.add_argument("-o",
                        "--output",
                        action='store',
                        nargs=1,
                        type=str,
                        help='The name of the output file',
                        default=None,
                        dest="output")



    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    print(args)
    #raise KeyboardInterrupt

    #if not len(args)==2:
    #    multi_epoch_diff(args.single_epoch_file[0], args.all_epoch_file[0])
    #else:
    multi_epoch_diff(args.single_epoch_file[0], args.all_epoch_file[0], args.output[0])
