'''
This module is a central repository of the values that the various contour
levels should have. Each script will check first in the directory from which
they are run for a personalized version of these files or default to the one
in pywrf/wrf/plot_wrfout/plot_wrfout_config.py
'''

import numpy as n

# I think only one style should be adopted for all the plots hence
plot_colorbar = True
plot_contour_labels = False
plot_wind_vectors = True
monochrome = False

frame_width = {}
frame_width['d01'] = 5.0
frame_width['d02'] = 5.0
frame_width['d03'] = 0.5

quiv_skip = {}
quiv_skip['d01'] = 5
quiv_skip['d02'] = 12
quiv_skip['d03'] = 10

meridians_delta = {}
meridians_delta['d01'] = 15
meridians_delta['d02'] = 15
meridians_delta['d03'] = 3

parallels_delta = {}
parallels_delta['d01'] = 15
parallels_delta['d02'] = 15
parallels_delta['d03'] = 3

quiverkey_length = {}
quiverkey_length['d01'] = 10
quiverkey_length['d02'] = 10
quiverkey_length['d03'] = 10

mslp_cntr_lvl = {}
# this looked pretty with the gist_ncar colormap
# over Australia at 36km resolution
# and grid corners (following wrfout metadata notation)
# LAT_LL_T = -50.81659
# LON_LL_T = 100.8601
# LAT_UR_T = -13.36167
# LON_UR_T = 165.0748
mslp_cntr_lvl['d01'] = n.arange(976,1041,4)
mslp_cntr_lvl['d02'] = n.arange(976,1041,4)
mslp_cntr_lvl['d03'] = n.arange(1004,1028,1)


sfc_mixing_ratio_cntr_lvl = {}
sfc_mixing_ratio_cntr_lvl['d01'] = n.arange(0,32,2)
sfc_mixing_ratio_cntr_lvl['d02'] = n.arange(0,32,2)
sfc_mixing_ratio_cntr_lvl['d03'] = n.arange(0,28,1)
