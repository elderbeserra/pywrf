'''
Physical and other constants in one place to achieve consistency across my
scripts.
'''

import numpy as n

# International Standard Atmosphere -> ISA -> isa_
isa_levels = 7
# m
isa_lower_boundaries = n.array([0, 11, 20, 32, 47, 51, 71]) * 1e3
# NB 86km (geometric) corresponds to 84.852km (geopotential)
# I have not derived the above so I am not sure what formula for the height
# dependecy of g -> TODO verify 
isa_upper_boundaries = n.array([11, 20, 32, 47, 51, 71, 84.852]) * 1e3
# Pa = N m^-2 = kg m s^-2 m^-2 = kg m^-1 s^-2
isa_sea_level_pressure = 101325
# K
isa_sea_level_temperature = 288.15
# K/m 
isa_temperature_gradient = n.array([-6.5, 0, 1, 2.8, 0, -2.8, -2]) * 1e-3


# m^3 kg^-1 s^-2 and plus minus 0.001
gravitational_constant = 6.6742e-11
# m s^-2
gravitational_acceleration = 9.80665

# air
# common g mol-1 here kg mol-1
nitrogen_molecular_weight          = 28.01 * 1e-3
oxygen_molecular_weight            = 32 * 1e-3
carbon_dioxide_molecular_weight    = 44.01 * 1e-3
water_molecular_weight             = 18.02 * 1e-3
air_molecular_weight               = 28.97 * 1e-3
# J kg^-1 K^-1 (approximate)
dry_air_gas_constant = 287
# J K^-1 mol^-1
universal_gas_constant = 8.3143
# J Kg^-1 K^-1
air_specific_heat_constant_pressure = 1004
air_specific_heat_constant_volume = 717
