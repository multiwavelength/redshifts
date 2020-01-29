# Archive of redshifts

This package searches within online astronomy archives Vizier and NED to find all unique spectroscopic redshifts within a designated field. 

## Usage
The main functionality of the package can be used by calling the function:
'''
query.run_query(data_path, name, coords)
'''
where data_path is the location where to place the downloaded data, name is the identifier for the field of interest and
coords are Astropy coordinates that specify where to point the query.
 
For details on how to use the package, please see main file, [PIPELINE_online_redshift.py](PIPELINE_online_redshift.py)

## Features and limitations
- Flexible search within a radius of a given set of (RA, DEC) coordinates
- Uses column names and descriptions (including columns contents input, UCD in Vizier) to identify columns containing spectroscopic redshifts of velocities
- Weeds out photometric redshifts and duplicates and returns a unique list of "best" spectroscopic redshift measurements.

## Limitations
The package relies on the original authors correctly using the UCD and other column names. However, there are clear cases of misuse, where labels reserved for spectroscopic redshifts contained photometric redshift. To remedy this, the package contains a list of "banned" catalogs, which was compiled by hand by inspecting catalogues that obviously misused column names and labels. The package also contains a list of banned keywords, that attempt to remove column that might contain redshift information, but not actual spectroscopic redshift measurements for individual sources.

## Dependencies
- Tested in Python 3.8, but should work in Python 3.x
- Python packages: os, sys, glob, numpy, shutil, astropy, astroquery, collections; To install, it is easiest to use pip:
```
pip install package-name
```
- Other packages: Starlink [STILTS](http://www.star.bris.ac.uk/~mbt/stilts/); To install STILTS, either use the repository (tests on Linux Mint) or use instructions from website:
```
sudo apt install stilts
```
