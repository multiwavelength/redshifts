from astropy import units as u

# Name of the column used to set the origin of the measurement 
origin_name = 'Source catalog'

################################################################################
# # List of keywords to have and avoid in searching the Vizier database

# Keywords we want in the column description for a redshift column
wanted_keywords_redshift = ['Redshift', 'redshift']

# Catalogues that have been found to have issues.
# GLADE: it reports photometric redshifts as redshift, rendering the catalogues
#        useless
banned_catalogues = ['glade1', 'glade2']

# Keywords in the description of the columns that indicate the column is not 
# actually a redshift measurement, but something else
banned_keywords = ['cluster', 'Cluster', 'sigma', 'Sigma', 'Photometric', 'photometric', 
                   'origin', 'Origin', 'reference', 'Reference', 'citation', 'Citation', 
                   'number', 'Number', 'uncertainty', 'Uncertainty', 
                   'error', 'Error on', 'Error of', 'distance modulus', 
                   'quality code', 'Quality flag', 'nearest neighbor', 'of the neighbors',
                   'Type of', 'Source of', 'Source for', 'photo-z', 
                   'Display the Object type', 'of the center', 'of the group center',
                   'corrected to the CMB', 'qualifier', 'Reliability of', 
                   'Peak of', 'Flag', 'source flag', 'FWHM']

# Units in the columns which indicate the column is not a redshift measurement,
# but something else 
banned_units = ["Mpc", "Y:M:D", "mag"]

# List of characters in the column name that indicate the column does not 
# contain a spectroscopic redshift measurement
banned_names = ["ph", "f_", "n_"]

# If the description contains these keywords, select this column over any other
# one
hard_selection = ["spectroscopic", "Spectroscopic"]

# Velocity unit for column mislabelled as redshift, when they actually contain
# velocities
vel_unit = u.km/u.s
