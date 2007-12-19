'''
Useful met formulas
'''
import numpy as n
import constants as c

def goff_gratch(temperature):
    '''
    computes saturation water pressure from temperature (in Kelvin)
    ref Smithsonian Tables 1984, after Goff and Gratch 1946
    usage:
    goff_gratch(temperature)
    '''
    exponent = \
      - 7.90298 * (373.16/temperature - 1) \
      + 5.02808 * n.log10(373.16 / temperature) \
      - 1.3816e-7 * (10**(11.344 * (1-temperature/373.16))  -1) \
      + 8.1328e-3 * (10**(-3.49149 * (373.16/temperature-1))  -1) \
      + n.log10(1013.246) 
    return 10**exponent

def water_vapour_pressure_to_mix_ratio(pressure, water_vapour_pressure):
    '''
    computes the mixing ratio given the total pressure and the water vapour
    partial pressure
    usage:
    water_vapour_pressure_to_mix_ratio(pressure, water_vapour_pressure)
    '''
    coeff = c.water_molecular_weight / c.air_molecular_weight
    dry_air_pressure = pressure - water_vapour_pressure
    return coeff * water_vapour_pressure / dry_air_pressure
