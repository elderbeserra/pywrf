from viz.utils import peek
from viz.utils import cardinal_2_month
import pylab as p
import matplotlib.numerix.ma as ma
import numpy as n
from string import zfill

def plot_soilw_from_gfs(file_name, cntr_lvl=None):
    file, vars = peek(file_name, show_vars=False)
    lon = vars['lon_3'].get_value()
    lat = vars['lat_3'].get_value()
    soilw_var = vars['SOILW_3_DBLY_10']
    soilw = soilw_var.get_value()
    soilw_m = ma.masked_where(soilw == soilw_var._FillValue, soilw)
    for lvl_idx in range(len(soilw)):
        p.figure()
        if cntr_lvl is not None:
           p.contourf(lon,lat,soilw_m[lvl_idx],cntr_lvl)
        else:
           p.contourf(lon,lat,soilw_m[lvl_idx])
        p.colorbar()
        # flip the y-axis
        p.ylim((lat[-1], lat[0]))
        p.ylabel('degrees North')
        p.xlabel('degrees East')
        lvl_bounds = vars['lv_DBLY5_l' + str(lvl_idx)]
        valid_when = soilw_var.initial_time
        valid_when = valid_when[3:5] + ' ' \
          + cardinal_2_month(int(valid_when[:2])) + ' ' + valid_when[6:10] \
          + ' ' + valid_when[12:17] + ' UTC'
        title_string = 'Volumetric soil moisture (fraction) valid at ' \
          + '\n' + valid_when + ' ' \
          + str(lvl_bounds[0]) + '-' + str(lvl_bounds[1]) + ' cm from GFS'
        p.title(title_string)
    return 

def plot_soilt_from_gfs(file_name, cntr_lvl=None):
    file, vars = peek(file_name, show_vars=False)
    lon = vars['lon_3'].get_value()
    lat = vars['lat_3'].get_value()
    soilt_var = vars['TMP_3_DBLY_10']
    soilt = soilt_var.get_value()
    soilt_m = ma.masked_where(soilt == soilt_var._FillValue, soilt)
    for lvl_idx in range(len(soilt)):
        p.figure()
        if cntr_lvl is not None:
           p.contourf(lon,lat,soilt_m[lvl_idx],cntr_lvl)
        else:
           p.contourf(lon,lat,soilt_m[lvl_idx])
        p.colorbar()
        # flip the y-axis
        p.ylim((lat[-1], lat[0]))
        p.ylabel('degrees North')
        p.xlabel('degrees East')
        lvl_bounds = vars['lv_DBLY5_l' + str(lvl_idx)]
        valid_when = soilt_var.initial_time
        valid_when = valid_when[3:5] + ' ' \
          + cardinal_2_month(int(valid_when[:2])) + ' ' + valid_when[6:10] \
          + ' ' + valid_when[12:17] + ' UTC'
        title_string = 'Soil temperature (K) valid at ' \
          + '\n' + valid_when + ' ' \
          + str(lvl_bounds[0]) + '-' + str(lvl_bounds[1]) + ' cm from GFS'
        p.title(title_string)
    return 

# plot_temp_from_gfs
# plot_wind_from_gfs
# plot_mr_from_gfs
# plot_ght_from_gfs
# plot_sfc_temp_from_gfs
# plot_sfc_wind_from_gfs
# plot_sfc_mr_from_gfs
# plot_mslp_from_gfs
# plot_sfc_press_from_gfs
# plot_land_from_gfs
# plot_terrain_from_gfs
# plot_skin_temp_from_gfs
