import yaml
from dataclasses import replace, asdict, is_dataclass, field
from typing import List

from pydantic.dataclasses import dataclass
from astropy import units as u

class Quantity(u.SpecificTypeQuantity):
    """
    Validation of the types of unit for each parameter, to ensure the right type
    is being given.
    """

    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, v):
        return cls(v)


class Angle(Quantity):
    _equivalent_unit = u.degree


@dataclass
class Setup:
    radius: Angle
    uncertainty: float
    banned_catalogs_redshift: List[str]
    banned_catalogs_velocity: List[str]


def read_config(config_file) -> Setup:
    """
    Read YAML configuration file into a class. Not all parameters have to be 
    set. It not set, a parameter will be set to the default value. The class has 
    defaults that the config file will override. The unit types will also be
    checked for correctness.
    Input:
        config_file: path to YAML parameter file
    Output:
        return a Constants dataclass instance with the defaults and config 
        overrides.
    """
    config = yaml.safe_load(open(config_file).read())
    return Setup(**config)
