import os
import shutil
import warnings
warnings.filterwarnings("ignore")
from astropy.io import fits
from astropy.table import Column, Table, QTable, vstack
import astropy.units as u
import astropy.coordinates as coord

import query as q
from astroquery.vizier import Vizier
import cross_match_utilities as cmu

for cluster in cluster_table[0:1]:
    coords = coord.SkyCoord(cluster['RA'], cluster['DEC'])
    
    gt = q.query_vizier(coords, 'redshift')
    gt.write('grand_Redshift.fits', format='fits', overwrite=True)

    gt = q.query_vizier(coords, 'velocity')
    gt.write('grand_Redshift_vel.fits', format='fits', overwrite=True)

    gt = q.query_NED(coords)
    gt.write('grand_Redshift_NED.fits', format='fits', overwrite=True)


import os
duration = 1  # seconds
freq = 440.  # Hz
os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq))
