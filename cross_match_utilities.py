import os, sys
import subprocess
import glob
import shutil

import numpy as np
from astropy.table import Column, Table
from collections import Counter
from collections import defaultdict

import setup_astroquery as sa


def add_origin(path1):
    """
    Add a column to zinfo.dat file with the name of the file, which thus
    indicate the origin of each of the redshift measurements
    Input:
        path1: zinfo.dat file
    Output:
        fits file with the extra column
    """
    table = rf.read_lof(path1)
    origin = Column([os.path.basename(path1)]*len(table), name='Origin', dtype='U', 
                          description='Origin of the redshift measurement')
    table.add_column(origin)
    outname = f'{path1.strip(".dat")}_origin.fits'
    table.write(outname, format='fits', overwrite=True)
    return outname


def concat_tables(list_of_tables, outname):
    """
    Concatenate using STILTS a number of fits tables.
    Input:
        list_of_tables: list of paths to each of the tables to be concatenated
        outname: preferred name for the outfile
    """
    no = len(list_of_tables)
    command = f'stilts tcatn nin={no} '
    for i, table in enumerate(list_of_tables):
        command = command + f'in{i+1}={table} '
    command = command + f'out={outname}'
    os.system(command)


def identify_duplicates(file1, outname, RA="RA", DEC="DEC", dist=1):
    """
    Within a single file,identify using STITLS whether there are any sources 
    with multiple redshift measurements
    Input: 
        file1: name of the file containing all the info on the targets. Must 
               contain RA and DEC under the names "RA" and "DEC"
        outname: preferred name for the output file. This will contain 2 extra
                 columns which describe whether the source has "duplicate" and
                 another which says how many duplicates exist in the group (this
                 will usually be 2, from 1 source with 2 measurements)
        RA, DEC: preferred names to be set for the RA and DEC column
        dist: distance in arcsec to be used for the crossmatch
    Output:
        return True if duplicates were found, False if no duplicates were found
        fits file with name outname which marks the duplicate sources   
    """
    command = (f'stilts tmatch1 matcher=sky values="{RA} {DEC}" params={dist} '+
               f'action=identify in={file1} out={outname}')
    try: 
        subprocess.check_output(command, shell=True, executable='/bin/zsh')
        return True
    except:
        shutil.copyfile(file1, outname)
        return False


def list_duplicates(seq):
    """
    Find groups of duplicated values in a sequence
    Input:
        seq: list of values
    Output:
        return a list of unique group of duplicate values 
    """
    tally = defaultdict(list)
    for i, item in enumerate(seq):
        try: 
            if item.mask==True: continue
        except:
            tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() 
                            if len(locs)>1)



def find_groups_redshift(file1, outfile, z=sa.z):
    """
    Search through the GroupID column of the fits table and list duplicates and
    their row id. The GroupID identifies pairs/triplets/groups of same sources
    which have multiple redshift estimations.
    Input:
        file1: master file with all the redshift for a cluster. May contain 
               duplicates from sources with multiple redshift measurements
    Output:
        fits file with unique sources 
    """

    table = Table.read(file1, hdu=1)
    ids = table['GroupID']
    indeces_remove = []

    # Iterate through the grouped sources
    for dup in sorted(list_duplicates(ids)):
        if dup[0]<0: continue # ignore negative ids which act as fillers for 
                              # unique sources with no matches
        # Find groups of sources
        grouped_sources = table[dup[1]]
    
        # Make a list of the length of the redshift as a proxy for the precision
        # of the redshift measurement
        significance = np.array([ len(format(source[z], 'f').rstrip('0')) 
                                                for source in grouped_sources ])
        # Remove the source with the most precision and add the rest to the list
        # of sources to be removed
        del dup[1][np.argmax(significance)]
        index_source_to_remove = dup[1]

        # Append all indeces to be removed from the list
        indeces_remove = indeces_remove+index_source_to_remove
    
    # Remove the lines corresponding to the sources chosen for removal through
    # the process described above
    table.remove_rows(np.array(indeces_remove))
    
    # Write out a new fits file containing only unique sources
    table.write(outfile, format='fits', overwrite=True)


def sum_SN(path_fits, line_list):
    """
    Calculate the summation of the signal-to-noise of all detected lines for a
    source
    Input:
        path_fits: path to the file containing the fits to the lines
        line_list: path to the file containing all the emission lines that were 
                   fit
    Output:
        sum of S/N for all detected lines
    """
    t = Table.read(path_fits, hdu=1)
    s = 0.
    for emline in line_list:
        line = emline['line']
        try:
            # Some lines might be completely out of the coverage area of the 
            # spectrum and can be ignored
            row = t[t.field('line')==line][0]
        except: 
            continue
        try: 
            if (np.abs(row['amplitude']/row['amplitude_err'])>const.SN_limit):
                s = s + np.abs(row['amplitude']/row['amplitude_err'])
        except:
            continue
    return s


def choose_unique_sources(path1, line_list):
    """
    Go through all of the *zinfo files, listing all the measured redshifts and 
    identify duplicate sources and produce a single list of unique sources,
    choosing the best S/N detections.
    Input:
        path1: path to the cluster in which all of the configurations are listed
        line_list: list of lines measured that are used for comparing the S/N
    Output:
        Unique list of sources with redshifts measured, written in a fits file 
    """
    dat_files = glob.glob(f'{path1}/*/*/*zinfo.dat')
    fits_files = []

    # Add origin/source of each target measurement to the table and save as 
    # a new fits file
    for dat_file in dat_files:
        fits_files.append(add_origin(dat_file))

    path_concat, path_ident, path_unique = u.naming_source_table(path1)

    # Concatenate all the lists of sources with redshifts into a single file
    concat_tables(fits_files, path_concat)
    
    # Find sources with multiple redshift measurements
    duplicates = identify_duplicates(path_concat, path_ident)
    # Write out a single fits file with unique sources, removing lower quality
    # duplicates
    if duplicates==True:
        find_groups(path_ident, path_unique, line_list)
    else:
        shutil.copyfile(path_ident, path_unique)
    

def cross_cat(file1, file2, outname, RA1, DEC1, RA2, DEC2, dist=1, join='all1'):
    """
    Cross match two catalogues, by default keeping all of the elements from the 
    first one as reference.
    Input: 
        file1: first file in the cross match, used as reference
        file2: second file in the cross match, only matches are selected
        outname: name of the resulting file to be written out
        RA, DEC: preferred names for the RA and DEC column in each table
        dist: distance in arcsec to be used for the crossmatch
    Output:
        return True if duplicates were found, False if no duplicates were found
        fits file with name outname which marks the duplicate sources   
    """
    command = (f'stilts tskymatch2 in1={file1} in2={file2} out={outname} ' +
               f'ra1={RA1} dec1={DEC1} ra2={RA2} dec2={DEC2} error={dist} '+
               f'join={join}')
    subprocess.check_output(command, shell=True, executable='/bin/zsh')