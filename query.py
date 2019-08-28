from astroquery.vizier import Vizier
from astroquery.ned import Ned
from astropy.table import Column, QTable, Table, vstack
from astropy import constants as const
from astropy import units as u

import setup_astroquery as sa
Ned.TIMEOUT = 100*60

def unwanted_catalogue(cat_name, banned_cat_list):
    """
    Remove unwanted catalogues, if their name matches the list of banned 
    catalogues
    Input: 
        cat_name: name of the catalogue
        banned_cat_list: list of catalogues that have been deemed banned
    Output:
        return True if the catalogue is within the banned list, False if not
    """
    if any(elem in cat_name for elem in banned_cat_list):
        return True
    else: 
        return False


def set_keywords(type):
    """
    Set the right collections of keywords wanted or banned within different 
    column names and descriptions. This helps select useful column and remove
    unwanted columns.
    Input:
        type: redshift of velocity
    Output:
        return the relevant collection of keywords (lists) located possibly in 
        another file
    """
    if type=='velocity':
        return sa.wanted_keywords_velocity, sa.banned_keywords_velocity, sa.banned_names_velocity
    if type=='redshift':
        return sa.wanted_keywords_redshift, sa.banned_keywords, sa.banned_names


def initial_check(type, desc, col_name):
    """
    Make initial checks that the columns fulfill a set of basic requirements. 
    Input:
        type: type of redshift search: velocity of redshift
        desc: description of the column
        col_name: column name
    Output:
        return True if the column fulfills all requirements; False if not
    """
    wk, bk, bn = set_keywords(type)
    if ( any(elem in desc for elem in wk) and 
        (not any(elem in desc for elem in bk)) and
        (not any(elem in col_name for elem in bn)) ):
        return True
    else:
        return False


def initial_select_columns(type, desc, col_name, unit=None):
    """
    Initial selection of column that contain the wanted keywords and none of the 
    unwanted keywords
    Input:
       desc: test description of the column
       col_name: column name
    """
    if initial_check(type, desc, col_name):
        if wanted_units(type, unit):
            return True
        else:
            return False
    else:
        return False


def wanted_units(type, unit):
    """
    Some columns have units which indicate they are a velocity measurement
    """
    if type=='velocity':
        if any(elem == unit for elem in sa.wanted_units_velocity):
            return True
        else: 
            return False
    if type=='redshift':
        return True


def banned_units(col_unit):
    """
    Some columns have units which indicate they are not a redshift measurement.
    Input:
        col_unit: a AStropy unit to test
    Output:
        If the unit in within that list return True; otherwise return False.
    """
    if ((col_unit is not None) and 
        (any(elem in col_unit.to_string() for elem in sa.banned_units)) ):
        return True
    else:
        return False


def vel2redshift(cat, col):
    """
    Sometimes a column is misidentified as a redshift, but it's really a 
    recessional velocity. Convert velocity to redshift and replace the column.
    Input: 
        cat: catalogue
        col: column name in string format
    Output:
        catalogue with velocity converted to redshift
    """
    if cat[col].unit==sa.vel_unit:
        desc = cat[col].info.description
        cat.replace_column(col, cat[col]/const.c.to(cat[col].unit))
        cat[col].info.description = desc
    return cat


def select_best_redshift(cat, col_selection):
    """
    Select a single column for the redshift measurement, for each catalogue. If 
    the table does not have any relevant column, return None. If there is one 
    relevant column, keep it. If presented with more columns that fill the 
    redshift requirement, if one of them is specifically designated as a 
    spectroscopic redshift, keep that one. Otherwise, just select the first 
    column in the list.
    Input:
        cat: parent catalogue
        col_selection: list of column names
    Return:
        None or a single column name
    """
    # If no columns were select, return None
    if len(col_selection)==0:
        return None

    # If exactly one column was selected, return its name
    if len(col_selection)==1:
        return col_selection[0]

    # If more than one column were previously selected, choose one.
    if len(col_selection)>1:
        for col in col_selection:
            if any(elem in cat[col].info.description for elem in sa.hard_selection):
                return col
        return col_selection[0]


def column_selection(type, cat):
    """
    Select the columns that could potentially contain redshift information
    Input:
        cat: catalogue in Astropy table format
    Output:
        selection of columns
    """
    col_selection = []
    for col in cat.colnames:
        desc = cat[col].info.description
        if banned_units(cat[col].unit):
            continue
        # Make first pass of selecting columns
        if initial_select_columns(type, desc, cat[col].info.name, cat[col].unit):
            # Ignore column that have units not consistent with redshift
            # measurements (like Mpc)
            # If column is velocity, convert to redshift
            cat = vel2redshift(cat, col)
            # Add valid column to the list of columns for further 
            # consideration
            col_selection.append(col)
    return col_selection


def set_unwanted_list(type):
    """
    Set which list of catalogues you want to ban. List can be in an outside file
    Input:
        type: velocity or redshift search
    Output:
        return list of catalogue names in string format
    """
    if type=='redshift':
        return sa.banned_catalogues
    if type=='velocity':
        return sa.banned_catalogues_velocity


def process_catalog(type, cat, RA, DEC, z, RAf, DECf):
    """
    Take an downloaded Vizier catalogue and extract RA, DEC and a redshift 
    column, if possible
    Input:
        type: velocity or redshift search
        cat: catalog as downloaded from Vizier with all columns in it
        RA, DEC: names of RA and DEC column in original catalog
        z, RAf, DECf: final names for the redshift, RA and DEC column
    Output:
        return None if catalog does not contain any useful redshift column or
        a table with 4 columns: RA, DEC, redshift and data origin
    """
    # Skip unwanted catalogues
    if unwanted_catalogue(cat.meta['name'], set_unwanted_list(type)): 
        return None
    
    # Make a list of potential column that contain redshift information
    col_selection = column_selection(type, cat)
    final_z_col = select_best_redshift(cat, col_selection)
    
    # If no relevant redshift column is present, skip the catalog
    if final_z_col == None: 
        return None
    # If all values are masked, skip the catalog
    if all(cat[final_z_col].mask): 
        return None

    # Homogenize column names
    cat.rename_column(final_z_col, z)
    
    # Select only relevant columns: RA, DEC and redshift
    final_cat = cat[RA, DEC, z][~cat[z].mask]

    # Rename the coord columns with chosen names
    final_cat.rename_column(RA, RAf)
    final_cat.rename_column(DEC, DECf)
    
    # Add Vizier catalog name to the table for future reference
    final_cat.add_column(Column([cat.meta['name']]*len(final_cat)), 
                         name=sa.origin_name)
    
    # Add to master list of tables
    return final_cat


def remove_potential_photoz(table, z_col='z_spec'):
    """
    Some redshifts in Vizier might still be photometric. Check whether the 
    values of the redshift has little precision and use this as a proxy for 
    labelling sources as photometric. Remove all rows consistent with 
    photometric measurements.
    Input:
        table: table that contains redshift measurements
        z_col: optional name of the redshift column 
    Return:
        table without potential photometric redshift
    """
    # Add number of rows to a list
    rows_to_remove = []
    # Iterate the table
    for i, z  in enumerate(table[z_col]):
        # A string of 9's or 0's usually comes from python not being able to 
        # represent float accurately, so those measurements probably have a lot
        # let precision than they might look like. Also if the string 
        # representation of the redshift measurement is short, it means the
        # measurement does not have a lot of precision; here I remove anything
        # with less than 2 significant digits.
        if ('999999' in str(z)) | ('000000' in str(z)) | (len(str(z))<=4):
            rows_to_remove.append(i)
    
    table.remove_rows(rows_to_remove)
    return table

def query_vizier(name, type, radius=sa.radius, RA='_RAJ2000', DEC='_DEJ2000', 
                 RAf=sa.RA, DECf=sa.DEC, z=sa.z):
    """
    Use astroquery to query the Vizier catalogue database
    Input:
        name: name of the source to query region for
        radius: optional, sets the radius to patrol around maine source
        RA, DEC: optional, coordinates of RA and DEC columns; defaults to vizier
                 names which point to RA and DEC homogenized to deg and J2000
        z: name to use for homogenized redshift column
    Return:
        Final table containing 4 columns, RA, DEC, redshift and origin of 
        redshift measurement, compiled from all data available of Vizier
    """

    # Setup the column and keywords used for the selection; calculate 
    # homogenised RA and DEC and return unlimited rows
    if type=='velocity':
        v = Vizier(columns=["**", RA, DEC], ucd="spect.dopplerVeloc*|phys.veloc*", 
                   row_limit=-1)
    if type=='redshift':
        v = Vizier(columns=["**", RA, DEC], ucd="src.redshift*", row_limit=-1)

    # Query a region using source name
    cat_list = v.query_region(name, radius=radius)

    # Initialize a list of table to be appended; These will be table that have
    # a redshift measurement
    table_list = []

    # Find whether the catalogue contains relevant information and save the 
    # column with the relevant data
    for cat in cat_list:
        final_cat = process_catalog(type, cat, RA, DEC, z, RAf, DECf)
        if final_cat is not None:
            table_list.append(final_cat)

    # Stack all the final catalogues with the relevant data
    grand_table = vstack(table_list)

    # Remove potential photoz's by removing measurements with little precision
    grand_table = remove_potential_photoz(grand_table)

    # Return results of the vizier search
    return grand_table


def fix_coord_units(cat, RA, DEC):
    """
    NED produces RA and DEC column with unrecognized astropy units. Fix them to 
    actual degrees.
    Input:
        cat: catalogue in question in Astropy table format
        RA, DEC: string names of the columns in question
    Return:
        catalogue with correct units
    """
    cat[RA].info.unit = u.deg
    cat[DEC].info.unit = u.deg
    return cat


def filter_ned_cat(cat, RA, DEC):
    """
    Filter the NED catalogue to include only sources that have redshift 
    measurements, remove photometric redshifts labelled as such and remove 
    groups of sources like clusters. Fix units of RA and DEC.
    Input:
        cat: catalogue in question in Astropy table format
        RA, DEC: string names of the columns for which units need to be fixed
    Return:
        catalogue with correct units, with rows containing redshift and without
        clusters
    """
    # Fix units
    cat = fix_coord_units(cat, RA, DEC)
    # Select only lines that have a redshift measurement
    cat = cat[~cat['Redshift'].mask]
    # Remove Galaxy Clusters
    cat = cat[~(cat['Type']==b'GClstr')]
    # Remove redshift values labelled as photometric
    cat = cat[~(cat['Redshift Flag']==b'PHOT')]
    return cat


def redshift_type(line, RA, DEC, un=sa.uncertainty):
    """
    Determine which type of redshift a source has associated with it. We do this
    by doing a targeted search on each source and checking its redshift 
    measurements. If there is no such table, or if the error on the redshift
    measurement is large, then the measurement is either missing or there are
    very large errors (ie the redshift is photometric). Return the line if the 
    values are good
    Input:
        line: line in original NED batch search
    Return
        return line of the table if it contains good spectroscopic redshift; 
        else, return None
    """
    try: 
        result_table = Ned.get_table(line['Object Name'], table='redshifts')
    except:
        return None
    if any(result_table['Published Redshift Uncertainty']<un):
        return line[RA, DEC, 'Redshift']
    else:
        return None


def query_NED(name, radius=sa.radius, RA='RA', DEC='DEC', z='z_spec', 
                                      RAf=sa.RA, DECf=sa.DEC):
    """
    Use astroquery to query the NED database
    Input:
        name: name of the source to query region for
        radius: optional, sets the radius to patrol around maine source
        RA, DEC: optional, coordinates of RA and DEC columns; defaults to vizier
                 names which point to RA and DEC homogenized to deg and J2000
        RAf, DECf, z: name to use for homogenized redshift and coord columns
    Return:
        Final table containing 4 columns, RA, DEC, redshift and origin of 
        redshift measurement, compiled from all data available of Vizier
    """
    # Query NED for the region around source within radius
    cat = Ned.query_region(name, radius=radius)
    
    # Filter the catalog, to remove useless rows
    filtered_cat = filter_ned_cat(cat, RA, DEC)
    
    # Initialize a list of table to be appended; These will be table that have
    # a redshift measurement
    row_list = []
    
    # Do another NED targeted search on each of the targets to check what type
    # of redshift it has associated
    for line in filtered_cat:
        if redshift_type(line, RA, DEC) is not None:
            row_list.append(line[RA, DEC, 'Redshift'])
    
    # Vstack all the lines into a single catalogue
    final_cat = vstack(row_list)
    
    # Rename columns to match general choice
    final_cat.rename_column('Redshift', z)
    final_cat.rename_column(RA, RAf)
    final_cat.rename_column(DEC, DECf)
    
    # Add origin column as NED
    final_cat.add_column(Column(['NED']*len(final_cat)), name=sa.origin_name)

    return final_cat

def query_redshift(target):
    """
    Perform an astroquery search of NED and Vizier for spectroscopic redshift 
    measurements. 
    Input: 
        target: either source name in string format or Astropy coordinate object
    Return:
        stacked table with all redshift measurements. Will most likely contain 
        duplicated sources
    """
    # Query Vizier for redshift columns
    tab_redshift = query_vizier(target, 'redshift')
    # Query Vizier for velocity columns
    tab_velocity = query_vizier(target, 'velocity')
    # Query NED for redshifts
    tab_NED = query_NED(target)

    return vstack([tab_redshift, tab_velocity, tab_NED])
