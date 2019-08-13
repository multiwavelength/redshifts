import warnings
warnings.filterwarnings("ignore")

import query as q


gt = q.query_vizier('Abell 115')
print(gt)
gt.write('grand_Redshift.fits', format='fits', overwrite=True)

import os
duration = 1  # seconds
freq = 440.  # Hz
os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq))