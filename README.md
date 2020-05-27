# **redshifts**

## What is **redshifts**

**redshifts** is a Python package which collects all unique spectroscopic redshifts from online databases such as Vizier and NED. 

## Features
- Flexible search within a radius of a given set of (RA, DEC) coordinates
- Uses column names and descriptions (including columns contents input, UCD in Vizier) to identify columns containing spectroscopic redshifts or velocities
- Weeds out photometric redshifts and duplicates and returns a unique list of "best" spectroscopic redshift measurements.

## Usage

To run **redshifts**, the following are needed:
- A configuration file, in YAML format, to specify things like the search radius, column names and banned catalogs
Details on the input file can be found further down. An example configuration file is contained within the git repository.


The output of **redshifts** include:
- A fits table with unique spectroscopic redshifts

**redshifts** can be used stand-alone from the terminal or from Python.
From the terminal, **redshifts** has a number of optional command line arguments. For details type:
```
redshifts --help
```
 

The main functionality of the package can be used by calling the function:
```
from redshifts.main import redshifts

redshifts(path, name, RA, DEC, config_file)
```
where path is the location where to place the downloaded data, name is the identifier for the field of interest, RA and DEC are coordinates that specify where to point the query and config_file is the YAML configuration file.


### Configuration file

The configuration file enables the user to customize the search.

A minimal working configuration example:
```yaml
radius: 0.01 degree
uncertainty: 0.002
banned_catalogs_redshift:
  - glade1
  - glade2
banned_catalogs_velocity:
  - VIII/7A/catalog
  - VIII/11/catalog
```

The search will be performed from the specified RA & DEC position out to the radius from the configuration file. The uncertainty is used to evaluate whether redshifts could be photometric instead of spectroscopic. The banned catalogs encompass any Vizier catalogs one does not want included in the search, for example because they were found to mix spectroscopic and photometric redshifts in one column.

## Limitations
The package relies on the original authors correctly using the UCD and other column names. However, there are clear cases of misuse, where labels reserved for spectroscopic redshifts contained photometric redshift. To remedy this, the package contains a list of "banned" catalogs, which was compiled by hand by inspecting catalogues that obviously misused column names and labels. The package also contains a list of banned keywords, that attempt to remove column that might contain redshift information, but not actual spectroscopic redshift measurements for individual sources.


## Installation requirements

- Python 3.8
- pip ^20.0
- Other packages: Starlink [STILTS](http://www.star.bris.ac.uk/~mbt/stilts/); To install STILTS, either use the repository (tested on Linux Mint) or use instructions from website:
```
sudo apt install stilts
```

## How to install **redshifts**
```
pip install git+https://github.com/multiwavelength/redshifts
```

