import os
import warnings

warnings.filterwarnings("ignore")
import click

import query as q
import constants as c


def redshifts(path, name, RA, DEC, config_path):
    # Read in configuration
    config = c.read_config(config_path)

    # Run the query
    q.run_query(path, name, RA, DEC, config)


@click.command()
@click.option("--path", default=".", help="Path for the data")
@click.option(
    "--config", default="redshiftconfig.yaml", help="Configuration file in YAML format."
)
@click.option("--name", default="", help="Target name.")
@click.option(
    "--RA", default="", help="Right ascension, with units: e.g. 150d, 150deg, 12h"
)
@click.option("--DEC", default="", help="Declination, with units: e.g. 30d, 30deg")
def main(path, config, name, ra, dec):
    # A fits file containing all the sources we want to download redshifts for
    redshifts(path, name, ra, dec, config)


if __name__ == "__main__":
    main()
