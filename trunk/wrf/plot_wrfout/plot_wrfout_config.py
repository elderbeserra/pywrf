'''
This module is a central repository of the values that the various contour
levels should have. Each script will check first in the directory from which
they are run for a personalized version of these files or default to the one
in pywrf/wrf/plot_wrfout/plot_wrfout_config.py
'''

import numpy as n
mslp_cntr_lvl = {}
# this looked pretty with the gist_ncar colormap
# over Australia at 36km resolution
# and grid corners (following wrfout metadata notation)
# LAT_LL_T = -50.81659
# LON_LL_T = 100.8601
# LAT_UR_T = -13.36167
# LON_UR_T = 165.0748
mslp_cntr_lvl['d01'] = n.arange(976,1041,4)
