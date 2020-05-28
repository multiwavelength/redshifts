import os
from io import BytesIO, StringIO
import shutil
import subprocess
from collections import defaultdict
import operator
from functools import reduce

from colorama import Fore
import numpy as np
from astroquery.vizier import Vizier
from astroquery.ned import Ned
from astropy.table import Column, QTable, Table, vstack
from astropy.io.votable import parse
from astropy import constants as const
from astropy import units as u
import astropy.coordinates as coord

timeout = 100 * 60
Ned.TIMEOUT = timeout

VELOCITY_SRC = "spect.dopplerVeloc*|phys.veloc*"
VELOCITY_KEYS = [
    "VELOC_BARYCENTER",
    "VELOC_HC",
    "VELOC_CMB",
    "VELOC_LG",
    "VELOC_LSR",
    "VELOC_GC",
]
REDSHIFT_SRC = "src.redshift*"
REDSHIFT_KEYS = ["REDSHIFT_HC"]
NED_FLAGS = [
    b"::",
    b"?",
    b"CONT",
    b"EST",
    b"FoF",
    b"LUM",
    b"MFA",
    b"MOD",
    b"PHOT",
    b"PEAK",
    b"PRED",
    b"SED",
    b"TENT",
    b"TOMO",
]
NED_TYPES = [b"QGroup", b"GClstr", b"GGroup", b"GPair", b"GTrpl", b"Other", b"PofG"]

HARD_SELECTION = ["spectroscopic", "Spectroscopic"]

BANNED_KEYWORDS = ["cluster", "Cluster", "Photometric", "photometric"]


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
    try:
        cat[col].unit.to(u.m / u.s)
        cat.replace_column(col, cat[col] / const.c.to(cat[col].unit))
        cat[col].info.description = (
            cat[col].info.description + ", converted to redshift"
        )
        cat[col].unit = u.dimensionless_unscaled
    except:
        pass
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
    # If no columns were selected, return None
    if len(col_selection) == 0:
        return None

    # If exactly one column was selected, return its name
    if len(col_selection) == 1:
        return col_selection[0]

    # If more than one column were previously selected, choose one.
    if len(col_selection) > 1:
        for col in col_selection:
            if any(elem in cat[col].info.description for elem in HARD_SELECTION):
                return col
        return col_selection[0]


def column_selection(type1, cat):
    """
    Select the columns that could potentially contain redshift information
    Input:
        type1: redshift of velocity
        cat: catalogue in Astropy table format
    Output:
        selection of columns
    """
    col_selection = []
    for col in cat.colnames:
        if col == "_RAJ2000":
            continue
        if col == "_DEJ2000":
            continue
        desc = cat[col].info.description
        f = any([(ban in desc) for ban in BANNED_KEYWORDS])
        if f is False:
            col_selection.append(col)
    return col_selection


def set_unwanted_list(type1, config):
    """
    Set which list of catalogues you want to ban. List can be in an outside file
    Input:
        type: velocity or redshift search
        c: configuration
    Output:
        return list of catalogue names in string format
    """
    if type1 == "redshift":
        return config.banned_catalogs_redshift
    if type1 == "velocity":
        return config.banned_catalogs_velocity


def process_catalog(
    type1, cat, config, RA, DEC, z="Redshift", RAf="RA", DECf="DEC", origin="Origin"
):
    """
    Take an downloaded Vizier catalogue and extract RA, DEC and a redshift 
    column, if possible
    Input:
        type1: velocity or redshift search
        cat: catalog as downloaded from Vizier with all columns in it
        config: configuration
        RA, DEC: names of RA and DEC column in original catalog
    Output:
        return None if catalog does not contain any useful redshift column or
        a table with 4 columns: RA, DEC, redshift and data origin
    """
    # Skip unwanted catalogues
    if unwanted_catalogue(cat.meta["name"], set_unwanted_list(type1, config)):
        return None

    # Make a list of potential column that contain redshift information
    col_selection = column_selection(type1, cat)
    final_z_col = select_best_redshift(cat, col_selection)

    # If no relevant redshift column is present, skip the catalog
    if final_z_col == None:
        return None

    # Skip weird column/tables with weird units/types for
    if cat[final_z_col].dtype not in [np.float32, np.float64]:
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
    final_cat.add_column(Column([cat.meta["name"]] * len(final_cat)), name=origin)

    # Add to master list of tables
    return final_cat


def remove_potential_photoz(table, z_col):
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
    for i, z in enumerate(table[z_col]):
        # A string of 9's or 0's usually comes from python not being able to
        # represent float accurately, so those measurements probably have a lot
        # less precision than they might look like. Also if the string
        # representation of the redshift measurement is short, it means the
        # measurement does not have a lot of precision; here I remove anything
        # with less than 2 significant digits.
        if ("999999" in str(z)) | ("000000" in str(z)):
            rows_to_remove.append(i)

    table.remove_rows(rows_to_remove)
    return table


def prelim_selection(cat_vot, type1, RA, DEC, KEYS):
    """
    Go through all the catalogs found online and decide whether to keep the
    catalog, and if yes, which columns.
    """
    cat_list = []
    for resource in cat_vot.resources:
        for table in resource.tables:
            cols = []
            for f in table.fields:
                # For velocities, check whether the units are the right type
                if type1 == "velocity":
                    if f.ucd in KEYS:
                        try:
                            f.unit.to(u.m / u.s)
                            cols.append(f.name)
                        except:
                            pass
                # For redshifts, check whether the precision is high enough
                if type1 == "redshift":
                    try:
                        if (
                            f.precision is not None
                            and (f.ucd in KEYS)
                            and (int(float(f.precision)) > 3)
                        ):
                            cols.append(f.name)
                    except:
                        pass
            if cols != []:
                tab1 = table.to_table(use_names_over_ids=True)[[RA, DEC] + cols]
                tab1 = tab1[reduce(operator.or_, [~tab1[col].mask for col in cols])]
                if type1 == "velocity":
                    for col in cols:
                        tab1 = vel2redshift(tab1, col)
                if len(tab1) > 0:
                    cat_list.append(tab1)
    return cat_list


def query_vizier(name, type1, config, RA="_RAJ2000", DEC="_DEJ2000", z="Redshift"):
    """
    Use astroquery to query the Vizier catalogue database
    Input:
        name: name of the source to query region for
        type1: redshift of velocity query
        config: configuration
        RA, DEC: optional, coordinates of RA and DEC columns; defaults to vizier
                 names which point to RA and DEC homogenized to deg and J2000
    Return:
        Final table containing 4 columns, RA, DEC, redshift and origin of 
        redshift measurement, compiled from all data available of Vizier
    """

    # Setup the column and keywords used for the selection; calculate
    # homogenized RA and DEC and return unlimited rows
    if type1 == "velocity":
        UCD = VELOCITY_SRC
        KEYS = VELOCITY_KEYS
    if type1 == "redshift":
        UCD = REDSHIFT_SRC
        KEYS = REDSHIFT_KEYS

    v = Vizier(columns=["**", RA, DEC], ucd=UCD, row_limit=-1, timeout=timeout,)

    # Query a region using source name, return a XML response and process
    # through Astropy as VOTABLE
    cat_vot = parse(
        BytesIO(v.query_region_async(name, radius=config.radius).text.encode()),
        pedantic=False,
        invalid="mask",
    )

    # Make a preliminary selection of columns to keep only RA, DEC and the
    # possible redshift columns
    cat_list = prelim_selection(cat_vot, type1, RA, DEC, KEYS)

    # Initialize a list of table to be appended; These will be tables that have
    # a redshift measurement
    table_list = []
    # Find whether the catalogue contains relevant information and save the
    # column with the relevant data
    for cat in cat_list:
        final_cat = process_catalog(type1, cat, config, RA, DEC)
        if final_cat is not None:
            table_list.append(final_cat)

    # Stack all the final catalogues with the relevant data
    if table_list:
        grand_table = vstack(table_list)

        # Remove potential photoz's by removing measurements with little precision
        grand_table = remove_potential_photoz(grand_table, z)

        # Return results of the vizier search
        return grand_table
    else:
        return None


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
    cat = cat[~cat["Redshift"].mask]

    # Remove Galaxy Clusters and Groups and redshifts labelled as photometric
    exclude_type = reduce(
        operator.and_, [(cat["Type"] != type_) for type_ in NED_TYPES]
    )
    exclude_flag = reduce(
        operator.and_, [(cat["Redshift Flag"] != flag) for flag in NED_FLAGS]
    )

    exclude = exclude_type & exclude_flag

    return cat[exclude]


def redshift_type(line, RA, DEC, un):
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
        result_table = Ned.get_table(line["Object Name"], table="redshifts")
    except:
        return None
    if any(result_table["Published Redshift Uncertainty"] < un):
        return line[RA, DEC, "Redshift"]
    else:
        return None


def query_NED(
    name,
    config,
    RA="RA",
    DEC="DEC",
    z="Redshift",
    RAf="RA",
    DECf="DEC",
    origin="Origin",
):
    """
    Use astroquery to query the NED database
    Input:
        name: name of the source to query region for
        config: configuration
        RA, DEC: optional, coordinates of RA and DEC columns; defaults to vizier
                 names which point to RA and DEC homogenized to deg and J2000    
    Return:
        Final table containing 4 columns, RA, DEC, redshift and origin of 
        redshift measurement, compiled from all data available of Vizier
    """
    # Query NED for the region around source within radius
    ned_result = Ned.query_region_async(name, radius=config.radius).text.encode()
    cat_vot = parse(BytesIO(ned_result), pedantic=False, invalid="mask")
    cat_vot = cat_vot.get_first_table().to_table(use_names_over_ids=True)

    # Filter the catalog, to remove useless rows
    filtered_cat = filter_ned_cat(cat_vot, RA, DEC)
    # Initialize a list of table to be appended; These will be table that have
    # a redshift measurement
    row_list = []

    # Do another NED targeted search on each of the targets to check what type
    # of redshift it has associated
    for line in filtered_cat:
        if redshift_type(line, RA, DEC, config.uncertainty) is not None:
            row_list.append(line[RA, DEC, "Redshift"])

    # Vstack all the lines into a single catalogue
    if not row_list:
        return None

    final_cat = vstack(row_list)

    # Rename columns to match general choice
    final_cat.rename_column("Redshift", z)
    final_cat.rename_column(RA, RAf)
    final_cat.rename_column(DEC, DECf)

    # Add origin column as NED
    final_cat.add_column(Column(["NED"] * len(final_cat)), name=origin)

    return final_cat


def query_redshift(target, path, name, config):
    """
    Perform an astroquery search of NED and Vizier for spectroscopic redshift 
    measurements. 
    Input: 
        target: either source name in string format or Astropy coordinate object
        path: path where to write the fits file with redshifts
        name: base name for the fits files containing the redshifts (can be the 
              target name)
        config: configuration
    Return:
        stacked table with all redshift measurements. Will most likely contain 
        duplicated sources
    """
    # Query NED for redshifts
    tab_NED = query_NED(target, config)
    if tab_NED is not None:
        tab_NED.meta["description"] = "NED"
        tab_NED.write(f"{path}/{name}/{name}_NED.fits", overwrite=True)
    print("Done querying NED...")
    # Query Vizier for redshift columns
    tab_redshift = query_vizier(target, "redshift", config)
    if tab_redshift is not None:
        tab_redshift.meta["description"] = "Vizier redshifts"
        tab_redshift.write(f"{path}/{name}/{name}_vizier_redshift.fits", overwrite=True)
    print("Done quering Vizier for redshifts...")
    # Query Vizier for velocity columns
    tab_velocity = query_vizier(target, "velocity", config)
    if tab_velocity is not None:
        tab_velocity.meta["description"] = "Vizier velocity"
        tab_velocity.write(f"{path}/{name}/{name}_vizier_velocity.fits", overwrite=True)
    print("Done querying Vizier for velocities...")

    # Add table to the list only if it not not empty
    cat_list = [t for t in [tab_redshift, tab_velocity, tab_NED] if t is not None]
    if cat_list:
        return vstack(cat_list)
    else:
        print(
            Fore.YELLOW + f"Warning: no spectroscopic redshifts found in search area."
        )
        exit()


def run_query(data_path, name, RA, DEC, config, z="Redshift"):
    """
    Query NED and Vizier for spectroscopic redshifts around coordinates of 
    choice
    Input: 
        data_path: location to place the downloaded data
        name: basename/identifier for field to query
        coords: coordinates in Astropy format
        config: configuration, including search radius, list of catalogues that 
                have been deemed banned
    Output:
        None; write out fits file with a unique list of redshifts around the
        target of interest
    """
    # Set the paths based on the where you want the data to be downloaded and
    # the identified for the sources/field the redshifts are downloaded for
    if not os.path.exists(f"{data_path}/{name}"):
        os.makedirs(f"{data_path}/{name}")
    path_concat = f"{data_path}/{name}/{name}_online_redshift.fits"
    path_ident = f'{path_concat.replace(".fits", "")}_ident.fits'
    path_unique = f'{path_concat.replace(".fits", "")}_ident_unique.fits'

    # Build coordinates
    coords = coord.SkyCoord(RA, DEC)

    # Perform the redshift query on Vizier and NED and write to fits file
    grand_table = query_redshift(coords, data_path, name, config)

    grand_table.meta["description"] = "Vizier and NED redshifts"
    grand_table.write(path_concat, format="fits", overwrite=True)

    # Identify duplicates and keep only the best redshift measurement
    duplicates = identify_duplicates(path_concat, path_ident, RA="RA", DEC="DEC")
    if duplicates == True:
        find_groups_redshift(path_ident, path_unique, z)
    else:
        shutil.copyfile(path_ident, path_unique)


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
    command = (
        f'stilts tmatch1 matcher=sky values="{RA} {DEC}" params={dist} '
        + f"action=identify in={file1} out={outname}"
    )
    try:
        subprocess.check_output(command, shell=True, executable="/bin/zsh")
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
            if item.mask == True:
                continue
        except:
            tally[item].append(i)
    return ((key, locs) for key, locs in tally.items() if len(locs) > 1)


def find_groups_redshift(file1, outfile, z):
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
    ids = table["GroupID"]
    indeces_remove = []

    # Iterate through the grouped sources
    for dup in sorted(list_duplicates(ids)):
        if dup[0] < 0:
            continue  # ignore negative ids which act as fillers for
        # unique sources with no matches
        # Find groups of sources
        grouped_sources = table[dup[1]]

        # Make a list of the length of the redshift as a proxy for the precision
        # of the redshift measurement
        significance = np.array(
            [len(format(source[z], "f").rstrip("0")) for source in grouped_sources]
        )
        # Remove the source with the most precision and add the rest to the list
        # of sources to be removed
        del dup[1][np.argmax(significance)]
        index_source_to_remove = dup[1]

        # Append all indeces to be removed from the list
        indeces_remove = indeces_remove + index_source_to_remove

    # Remove the lines corresponding to the sources chosen for removal through
    # the process described above
    table.remove_rows(np.array(indeces_remove))

    # Write out a new fits file containing only unique sources
    table.write(outfile, format="fits", overwrite=True)
