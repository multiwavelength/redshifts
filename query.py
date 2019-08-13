from astroquery.vizier import Vizier
from astroquery.ned import Ned
from astropy.table import Column, QTable, Table, vstack
from astropy import constants as const
from astropy import units as u

import setup_vizier as sv


def unwanted_catalogue(cat_name):
    """
    Remove unwanted catalogues, if their name matches the list of banned catalogues
    Input: 
        cat_name: name of the catalogue
    """
    if any(elem in cat_name for elem in sv.banned_catalogues):
        return True
    else: 
        return False


def initial_select_columns(desc, col_name):
    """
    Initial selection of column that contain the wanted keywords and none of the 
    unwanted keywords
    Input:
       desc: test description of the column
       col_name: column name
    """
    if ( any(elem in desc for elem in sv.wanted_keywords_redshift) and 
        (not any(elem in desc for elem in sv.banned_keywords)) and
        (not any(elem in col_name for elem in sv.banned_names)) ):
        return True
    else:
        return False


def banned_units(col_unit):
    """
    Some columns have units which indicate they are not a redshift measurement
    """
    if ((col_unit is not None) and 
        (any(elem in col_unit.to_string() for elem in sv.banned_units)) ):
        return True
    else:
        return False


def vel2redshift(cat, col):
    """
    Sometimes a column is misidentified as a redshift, but it's really a 
    recessional velocity. Convert to redshift and replace the column
    Input: 
        cat: catalogue
        col: column name in string format
    """
    if cat[col].unit==sv.vel_unit:
        cat.replace_column(col, cat[col]/const.c.to(cat[col].unit))
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
            if any(elem in cat[col].info.description for elem in sv.hard_selection):
                return col
        return col_selection[0]


def column_selection(cat):
    col_selection = []
    for col in cat.colnames:
        desc = cat[col].info.description
        # Make first pass of selecting columns
        if initial_select_columns(desc, cat[col].info.name):
            # Ignore column that have units not consistent with redshift
            # measurements (like Mpc)
            if banned_units(cat[col].unit):
                continue   
            else:
                # If column is velocity, convert to redshift
                cat = vel2redshift(cat, col)
                # Add valid column to the list of columns for further 
                # consideration
                col_selection.append(col)
    return col_selection


def process_catalog(cat, RA, DEC, z):
    """
    Take an downloaded Vizier catalogue and extract RA, DEC and a redshift 
    column, if possible
    Input:
        cat: catalog as downloaded from Vizier with all columns in it
    Output:
        return None if catalog does not contain any useful redshift column or
        a table with 4 columns: RA, DEC, redshift and data origin
    """
    # Skip unwanted catalogues
    if unwanted_catalogue(cat.meta['name']): 
        return None
    
    # Make a list of potential column that contain redshift information
    col_selection = column_selection(cat)
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
    
    # Add Vizier catalog name to the table for future reference
    final_cat.add_column(Column([cat.meta['name']]*len(final_cat)), 
                         name=sv.origin_name)
    
    # Add to master list of tables
    return final_cat


def query_vizier(name, radius=0.5*u.deg, RA='_RAJ2000', DEC='_DEJ2000', z='z_spec'):
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
    v = Vizier(columns=["**", RA, DEC], ucd="src.redshift*", row_limit=-1)

    # Query a region using source name
    cat_list = v.query_region(name, radius=radius)

    # Initialize a list of table to be appended; These will be table that have
    # a redshift measurement
    table_list = []

    # Find whether the catalogue contains relevant information and save the 
    # column with the relevant data
    for cat in cat_list:
        final_cat = process_catalog(cat, RA, DEC, z)
        if final_cat is not None:
            table_list.append(final_cat)

    # Stack all the final catalogues with the relevant data
    grand_table = vstack(table_list)
    
    # Return results of the vizier search
    return grand_table
