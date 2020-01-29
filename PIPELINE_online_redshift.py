import os
import warnings
warnings.filterwarnings("ignore")
from astropy.table import QTable
import astropy.coordinates as coord

import query as q
 
# Set the paths where the data will be downloaded
base_path = f'/home/andra/Desktop/Keep/astroquery'
data_path = f'{base_path}/clusters'
pipeline_path = f'{base_path}/pipeline'

# A fits file containing all the sources we want to download redshifts for
cluster_table = QTable.read(f'{data_path}/Cluster_properties.fits')

for cluster in cluster_table:
    # Set the base name/identifier for the source around which redshifts will
    # be searched for 
    name = cluster['Cluster']
    
    if not os.path.exists(f'{data_path}/{name}'):
        os.mkdir(f'{data_path}/{name}')
    
    # Set the central coordinates; in this case they are taken from the cluster
    # table as the RA and DEC of the cluster, but this can be set to anything
    coords = coord.SkyCoord(cluster['RA'], cluster['DEC'])
    
    # Set the paths based on the where you want the data to be downloaded and 
    # the identified for the sources/field the redshifts are downloaded for    
    path_concat = f'{data_path}/{name}/{name}_online_redshift.fits'
    path_ident = f'{path_concat.replace(".fits", "")}_ident.fits'
    path_unique = f'{path_concat.replace(".fits", "")}_ident_unique.fits'

    # Perform the redshift query on Vizier and NED and write to fits file
    grand_table = q.query_redshift(coords, data_path, name)
    grand_table.meta['description'] = u'Vizier and NED redshifts'
    grand_table.write(path_concat, format='fits', overwrite=True)

    # Identify duplicates and keep only the best redshift measurement
    duplicates = cmu.identify_duplicates(path_concat, path_ident, 
                                         RA='RA', DEC='DEC')
    if duplicates==True:
        cmu.find_groups_redshift(path_ident, path_unique)
    else:
        shutil.copyfile(path_ident, path_unique)
    
import os
duration = 1  # seconds
freq = 440.  # Hz
os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq))
