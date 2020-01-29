from astropy import units as u

# Search radius
radius = 0.7*u.deg

# Name of the column used to set the origin of the measurement 
origin_name = 'Origin'
# Final names for redshift, RA and DEC columns
z = 'Redshift'
RA = 'RA'
DEC = 'DEC'

################################################################################
# List of keywords to have and avoid in searching the Vizier database for 
# redshift

# Keywords we want in the column description for a redshift column
wanted_keywords_redshift = ['Redshift', 'redshift']

# Catalogues that have been found to have issues.
# GLADE: it reports photometric redshifts as redshift, rendering the catalogues
#        useless
banned_catalogues = ['glade1', 'glade2', 'IX/31/wgacat', 'VIII/29A/abelclus', 
                     'J/A+A/425/367/reflex50', 'J/A+A/425/367/reflex70', 
                     'VII/283/catalog', 'VII/258/vv10', 'IX/31/wgacat', 
                     'J/A+A/525/A157', 'J/ApJS/196/11', 'J/A+AS/103/349',
                     'J/ApJS/146/267', 'VII/259/6dfgs', 'VII/69/catalog',
                     'J/MNRAS/277/1477', 'J/MNRAS/478/1512', 'J/MNRAS/443/1231']

banned_catalogues_velocity = ['glade1', 'glade2', 'I/348/catalog', 'J/A+A/525/A157',
        'II/335/galex_ais', 'II/335/table4', 'V/139/sdss9', 'V/147/sdss12', 
        'V/148/morx', 'V/148/morxsupp', 'V/149/dr2', 'V/149/mstars2', 'V/153/dr4', 
        'V/153/astars4', 'V/153/mstars4', 'VII/177/table1', 'VII/190/zg_ori', 
        'VII/273/hmq', 'VII/279/dr12q', 'VII/279/dr12qbal', 'VII/283/catalog', 
        'VIII/3/catalog', 'VIII/7A/catalog', 'IX/10A/1rxs', 'IX/10A/1rxs_cor', 
        'IX/10A/cor_ned', 'IX/10A/cor_nvs', 'IX/10A/cor_ros', 'B/iram/30m', 
        'J/ApJ/448/521/table11', 'J/ApJ/511/65/table3', 'J/ApJ/554/L129/system', 
        'J/ApJ/790/158/table1', 'J/ApJ/798/73/table6', 'J/ApJ/817/112/highpm', 
        'J/ApJ/847/123/table1', 'J/ApJS/154/673/DIRBE', 'J/ApJS/235/11/table1', 
        'J/A+A/575/A48/table2', 'J/A+A/594/A116/nhi_hpx', 'J/A+A/596/A14/grlist_2', 
        'J/A+A/596/A14/grlist_s', 'J/A+A/596/A14/galist_2', 'J/A+A/596/A14/galist_s', 
        'J/A+A/619/L8/table2', 'J/AJ/149/171/table4', 'J/AJ/154/251/table8', 
        'J/MNRAS/416/2840/catalog', 'J/MNRAS/450/3675/scall', 
        'J/MNRAS/450/3675/scsingle', 'J/AZh/87/760/clusters', 'J/other/APh/26.282/tableb8',
        'VII/251/sbs', 'J/A+A/540/A106/dr8gr', 'J/ApJS/234/1/dwarfs', 
        'J/ApJ/687/78/table2', 'IX/10A/cor_ver', 'J/ApJS/172/599/table3', 
        'J/ApJ/819/63/table4', 'J/ApJ/680/169/table1', 'J/MNRAS/422/25/table1', 
        'V/149/astars2', 'J/A+A/596/A14/fpd_be_m', 'J/A+A/479/927/groups', 
        'J/A+A/532/A74/efigi', 'J/ApJS/194/45/catalog', 'J/AJ/149/171/table3', 
        'J/ApJ/614/91/table1', 'J/A+A/596/A14/fpd_be_s', 'J/AJ/126/1720/table1', 
        'J/ApJ/804/L15/data', 'J/A+A/596/A14/fpd_gr_m', 'J/A+A/602/A100/table1', 
        'J/A+A/479/927/galaxies', 'V/130/gcs3', 'J/AJ/142/122/table1', 
        'J/A+A/566/A1/galaxies', 'J/A+A/566/A1/v190_gr', 'J/ApJ/737/71/table1', 
        'J/A+A/596/A14/fpd_ga_s', 'J/A+A/596/A14/fpd_ga_p', 'J/ApJS/186/427/table2', 
        'J/MNRAS/410/860/table1', 'J/A+A/566/A1/v180_gr', 'J/MNRAS/377/787/lrgs', 
        'VII/256/table1', 'J/AJ/139/1808/table1', 'J/ApJS/165/1/table4', 
        'J/A+A/602/A100/table2', 'J/A+A/566/A1/v200_gr', 'J/ApJ/852/72/mockgals', 
        'J/A+A/566/A1/v195_gr', 'J/A+A/596/A14/fpd_ga_m', 'VII/259/spectra', 
        'J/A+A/596/A14/fpd_be_p', 'J/A+A/596/A14/fpd_gr_p', 'J/AJ/108/1/table1', 
        'J/MNRAS/446/3943/catalog', 'J/ApJ/778/188/table1', 'J/A+A/531/A165/param', 
        'J/A+A/566/A1/v205_gr', 'J/A+A/596/A14/fpd_gr_s', 'J/A+A/540/A106/dr8gal', 
        'J/MNRAS/439/611/catalog', 'J/ApJ/763/37/table2', 'J/MNRAS/455/2440/catalog', 
        'J/A+A/566/A1/groups', 'J/AJ/140/817/table2', 'J/AJ/142/188/table5', 
        'J/A+A/403/555/table2', 'J/A+A/562/A64/table3', 'J/ApJS/167/1/table4', 
        'J/MNRAS/462/2980/table1', 'J/ApJS/171/146/Stars', 'J/A+A/536/A43/table4s', 
        'J/ApJ/682/964/table4', 'J/A+A/578/A134/trujillo', 'J/ApJ/835/280/galaxies', 
        'J/A+A/377/801/table1', 'J/MNRAS/389/1074/tablea1', 'J/A+AS/106/303/table2', 
        'J/MNRAS/379/867/table1', 'J/MNRAS/392/135/clusters', 'J/AJ/146/151/table6', 
        'J/A+AS/111/237/table3', 'J/ApJ/729/22/table2', 'J/AJ/129/1063/table1', 
        'J/A+AS/139/545/table1', 'J/AJ/130/968/table2', 'J/ApJS/215/12/table1', 
        'J/A+A/412/657/table1', 'VII/269/dr9q', 'J/ApJS/225/23/table5', 
        'J/AJ/125/1817/table3a', 'J/AN/327/365/grGal', 'J/MNRAS/434/163/table2', 
        'J/ApJS/141/503/table1', 'VII/269/dr9qsup', 'J/A+A/575/A4/table4', 
        'J/A+A/551/A24/table3', 'J/ApJ/836/115/table3', 'VII/193/zbig', 
        'J/other/RAA/14.1135/table1', 'J/A+A/525/A90/tableb2', 
        'J/AJ/125/1817/table3d', 'J/AJ/125/1817/table2c', 'J/A+A/375/661/table2', 
        'J/MNRAS/413/1024/table1', 'J/MNRAS/372/1425/tablec2', 
        'J/MNRAS/443/1231/FPsample', 'J/MNRAS/366/144/table3', 
        'J/AJ/147/86/table6', 'J/ApJ/866/137/table2', 'J/AJ/149/27/table3', 
        'J/AJ/154/115/table2', 'J/A+A/548/A66/dr9qsup', 'VII/270/dr10qsup', 
        'J/A+A/348/897/tablea6', 'J/AJ/125/1817/table2d', 'J/A+A/588/A98/table3', 
        'J/A+A/551/A24/table1', 'J/A+A/578/A61/tableb2', 'J/MNRAS/449/1593/table1', 
        'J/A+A/525/A90/tableb1', 'J/ApJS/215/14/table4', 'J/PAZh/34/446/table2', 
        'J/ApJ/772/82/table1', 'I/325/icrf-u2', 'J/MNRAS/443/1231/table2', 
        'J/A+A/616/A114/tablea1', 'J/MNRAS/348/866/2dfsgpga', 'I/345/ssoobs', 
        'J/ApJS/225/23/table1', 'J/MNRAS/484/L8/stars', 'J/ApJS/233/3/table6', 
        'J/AJ/147/86/table4', 'J/MNRAS/412/1419/table2', 'J/MNRAS/348/866/2dfngpgr', 
        'J/ApJ/594/1/table8', 'J/AN/327/365/isoGal', 'J/ApJ/835/280/table3', 
        'J/ApJS/167/1/table3', 'J/A+AS/105/211/table2', 'J/AN/327/365/groups', 
        'J/AJ/156/18/table2', 'J/ApJ/711/284/table3', 'J/MNRAS/329/87/table1', 
        'J/A+A/423/683/table2', 'J/AJ/125/1817/table2b', 'J/ApJ/669/791/qsos', 
        'J/ApJS/225/23/table3', 'J/ApJ/763/32/table1', 'J/A+A/403/555/table1', 
        'J/A+A/436/443/table3', 'J/ApJ/705/1099/lenses', 'B/sb9/main', 
        'J/ApJ/849/20/eco', 'J/ApJS/225/23/table6', 'J/AJ/155/256/sample', 
        'J/ApJ/729/22/table3', 'J/PASJ/53/517/table3', 'J/ApJS/176/414/table1', 
        'J/AJ/125/1817/table2a', 'J/ApJ/744/177/table1', 'I/196/annex1', 
        'J/A+A/523/A91/table11', 'J/AJ/125/1817/table3c', 'J/AJ/153/257/table1c', 
        'J/ApJS/162/207/stars', 'J/A+A/523/A91/table10', 'J/A+A/545/A15/table1', 
        'J/MNRAS/410/190/table1', 'J/ApJS/210/5/table1', 'J/A+A/474/783/table2', 
        'J/ApJS/210/7/table1', 'J/A+A/520/A79/tablea1', 'J/ApJS/220/3/catalog', 
        'V/103/stars', 'J/MNRAS/348/866/2dfngpga', 'J/A+AS/135/511/table5', 
        'J/A+A/541/A40/table1', 'J/ApJ/743/171/2LAC', 'VII/270/dr10q', 
        'J/A+A/616/A10/tablea1a', 'IX/10A/cor_iras', 'J/ApJ/806/185/table1', 
        'J/MNRAS/372/1425/tableb1', 'J/AJ/147/86/table7', 'J/A+A/566/A1/v185_gr', 
        'J/MNRAS/443/1231/table9', 'J/MNRAS/449/2345/table1', 'J/ApJS/116/1/starcat', 
        'J/MNRAS/454/3962/table2', 'J/ApJS/221/32/table1', 'J/A+A/548/A66/dr9q', 
        'J/AJ/153/208/targets', 'J/A+A/566/A1/v210_gr', 'J/ApJ/690/670/table1', 
        'J/AJ/125/1817/table3b', 'J/MNRAS/448/2530/tablea1', 'J/MNRAS/348/866/2dfsgpgr', 
        'J/A+A/461/397/tablea1', 'J/ApJ/791/88/table1', 'J/ApJ/708/661/sn', 
        'V/125/obubvbet', 'J/A+A/564/A79/pm_ucac4', 'VIII/11/catalog', 
        'J/ApJ/698/819/table1', 'J/MNRAS/465/2120/table2', 'J/ApJS/196/11',
        'J/MNRAS/465/2120/tableb', 'J/ApJS/154/673/table5', 'J/A+AS/103/349',
        'J/ApJ/789/23/table2', 'J/ApJS/234/31/table11', 'V/125/obcat', 
        'J/A+A/309/749/tablea1', 'J/A+A/504/347/clusters', 'J/AJ/155/181/table1',
        'J/ApJS/146/267', 'J/MNRAS/488/590']

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
                   'Peak of', 'Flag', 'source flag', 'FWHM', 'Identifier']

# Units in the columns which indicate the column is not a redshift measurement,
# but something else 
banned_units = ["Mpc", "Y:M:D", "mag"]

# List of characters in the column name that indicate the column does not 
# contain a spectroscopic redshift measurement
banned_names = ["ph", "f_", "n_", "r_", "e_"]

# If the description contains these keywords, select this column over any other
# one
hard_selection = ["spectroscopic", "Spectroscopic"]

# Velocity unit for column mislabelled as redshift, when they actually contain
# velocities
vel_unit = u.km/u.s

################################################################################
# List of keywords to have and avoid in searching the Vizier database for 
# redshift

# Keywords we want in the column description for a redshift column
wanted_keywords_velocity = ['velocity', 'Velocity', 'redshift', 'Redshift']
wanted_units_velocity = [u.km/u.s]

# Keywords in the description of the columns that indicate the column is not 
# actually a redshift measurement, but something else
banned_keywords_velocity = ['proper motion', 'Unreliable', 'dispersion',
                            'expansion', 'Expansion', 'transverse', 'Transverse',
                            'Tangential', 'tangential', 'Rotational', 'rotational',
                            'peculiar', 'Peculiar', 'error of', 'Uncertainty of', 
                            'CIV troughs', 'BAL troughs', 'velocity width', 
                            'component of', 'radial velocity deviation',
                            'precision needed for', 'sigma', 'CMB frame', 
                            'CMB reference', 'microturbulence', 'Microturbulence',
                            'microturbulent', 'Microtubulent', 'apparent velocity', 
                            'velocity shift', 'Velocity shift', 'group velocity', 
                            'velocity of group', 'Group center velocity',
                            'component U', 'component V', 'component W',
                            'velocity U', 'velocity V', 'velocity W',
                            'U vel', 'V vel', 'W vel', 'Uvel', 'Vvel', 'Wvel',
                            'along x', 'along y', 'along z', 'Vx', 'Vy', 'Vz',
                            'VX', 'VY', 'VZ', 'center U', 'rotation V', 'Pole W',
                            'in local sheet reference frame', 'Interstellar',
                            'Local Standard of Rest', 'cylindrical reference frame',
                            'on GLON', 'on GLAT', 'radial velocity difference',
                            'Minimum velocity', 'Maximum velocity', 
                            'maximum velocity', 'minimum velocity',
                            'Radial velocity correction', 'velocity offset', 
                            'Velocity difference between', 'cluster velocity',
                            'cluster mean velocity']

# Units in the columns which indicate the column is not a redshift measurement,
# but something else 
banned_units_velocity = ["Mpc", "Y:M:D", "mag"]

# List of characters in the column name that indicate the column does not 
# contain a spectroscopic redshift measurement
banned_names_velocity =["e_"] #["ph", "f_", "n_"]

# If the description contains these keywords, select this column over any other
# one
hard_selection_velocity = ["spectroscopic", "Spectroscopic"]

################################################################################
# List of keywords to have and avoid in searching the NED database

# Redshift uncertainty for a redshift measurement to be considered spectroscopic
uncertainty = 0.002
