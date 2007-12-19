'''
Functions to calculate thermodynamic quantities for the International Standard
Atmosphere.
'''

from __future__ import division
from sys import path
path.append('/Users/val/projects/vapylib')
import numpy as n
import constants 

isa_levels = constants.isa_levels
lower_boundaries = constants.isa_lower_boundaries
upper_boundaries = constants.isa_upper_boundaries
sea_level_temperature = constants.isa_sea_level_temperature
sea_level_pressure = constants.isa_sea_level_pressure
temperature_gradient = constants.isa_temperature_gradient
g = constants.gravitational_acceleration
R = constants.dry_air_gas_constant

def which_level(geopotential_height):
    # basic check
    if geopotential_height < lower_boundaries[0] \
        or geopotential_height > upper_boundaries[-1]:
        message = 'z must be > ' + str(lower_boundaries[0]) + ' and < ' \
            + str(upper_boundaries[-1])
        raise 'BadGeopHgt', message
    for level_idx in range(isa_levels):
        if geopotential_height <= upper_boundaries[level_idx]:
            return level_idx
    
def temperature_calculator(geopotential_height):
    level = which_level(geopotential_height)
    temperature = sea_level_temperature
    for level_idx in range(level):
        temperature += temperature_gradient[level_idx] \
            * (upper_boundaries[level_idx] - lower_boundaries[level_idx])
    temperature += temperature_gradient[level] \
            * (geopotential_height - lower_boundaries[level])
    return temperature


def pressure_calculator(geopotential_height):
    # following Richard's notes (eqn 3.3)
    level = which_level(geopotential_height)
    integral = 0
    for level_idx in range(level):
        lower_boundary_temperature = \
            temperature_calculator(lower_boundaries[level_idx])
        if temperature_gradient[level_idx] != 0.:
            integral += n.log(  \
                (lower_boundary_temperature \
                + (upper_boundaries[level_idx] - lower_boundaries[level_idx]) \
                * temperature_gradient[level_idx]) / lower_boundary_temperature \
                ) / temperature_gradient[level_idx]
        else:
            integral += \
             (upper_boundaries[level_idx] -lower_boundaries[level_idx]) \
             / lower_boundary_temperature
    lower_boundary_temperature = \
        temperature_calculator(lower_boundaries[level])
    if temperature_gradient[level] != 0.:
        integral += n.log(  \
            (lower_boundary_temperature \
            + (geopotential_height - lower_boundaries[level]) \
            * temperature_gradient[level]) / lower_boundary_temperature \
            ) / temperature_gradient[level]
    else:
        integral += \
         (geopotential_height -lower_boundaries[level]) \
         / lower_boundary_temperature
    return sea_level_pressure * n.exp( -g / R * integral)

def find_height(pressure,tolerance=1):
    residue = 1e6
    top = upper_boundaries[-1]
    bottom = lower_boundaries[0]
    while abs(residue) > tolerance:
        midpoint = (top + bottom) / 2
        midpoint_pressure = pressure_calculator(midpoint)
        residue = pressure - midpoint_pressure
        #print midpoint, top, bottom, midpoint_pressure, residue, abs(residue), tolerance
        # assuming p decreases with height...
        if residue > 0:
            top = midpoint
        else:
            bottom = midpoint
    return int(midpoint)




