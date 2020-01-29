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

    # Run the query
    q.run_query(data_path, name, coords)