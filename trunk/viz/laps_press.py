from viz.utils import peek
from viz.utils import cardinal_2_month
from viz.utils import set_default_basemap
import pylab as p
import matplotlib.numerix.ma as ma
import numpy as n
from string import zfill

a_small_number = 1e-8

def plot_soilw(file_name, cntr_lvl=None):
    file, vars = peek(file_name, show_vars=False)
    lon = vars['lon'].get_value()
    lat = vars['lat'].get_value()
    soilw_var = vars['soil_mois']
    soilw = soilw_var.get_value()[0]
    mask = n.zeros(soilw.shape, dtype=int)
    for lvl_idx in range(len(soilw)):
        mask[lvl_idx] = vars['sea_lnd_ice_mask'].get_value()
    soilw_m = ma.masked_where(mask == 0 , soilw)
    lvl_bounds = vars['soil_lvl'].get_value()
    valid_date = str(vars['valid_date'].get_value()[0])
    valid_time = zfill(str(vars['valid_time'].get_value()[0]), 4)
    valid_when = valid_date[6:] + ' ' \
      + cardinal_2_month(int(valid_date[4:6])) + ' ' \
      + valid_date[0:4] \
      + ' ' + valid_time[:2] + ':' \
      + valid_time[2:] + ' UTC'
    m = set_default_basemap(lon,lat)
    # must plot using 2d lon and lat
    LON, LAT = p.meshgrid(lon,lat)
    #for lvl_idx in range(1):
    for lvl_idx in range(len(soilw)):
        p.figure()
        if cntr_lvl is not None:
            # let's not forget the scaling by the layer depth
            m.contourf(LON,LAT,soilw_m[lvl_idx]/lvl_bounds[lvl_idx], cntr_lvl)
        else:
            # let's not forget the scaling by the layer depth
            m.contourf(LON,LAT,soilw_m[lvl_idx]/lvl_bounds[lvl_idx])
        m.drawcoastlines()
        m.drawmeridians(n.array(n.arange(lon.min(), lon.max() + a_small_number, 15.)), labels=[1,0,0,1])
        m.drawparallels(n.array(n.arange(lat.min(), lat.max() + a_small_number, 15.)), labels=[1,0,0,1])
        p.colorbar(orientation='horizontal', shrink=0.7, fraction=0.02, pad=0.07, aspect=50)
        if lvl_idx == 0:
            title_string = 'Volumetric soil moisture (fraction) valid at ' \
              + '\n' + valid_when + ' ' \
              + '0' \
              + '-' + str(int(round(lvl_bounds[lvl_idx]*100))) + ' cm' \
              + ' from LAPS'
        else:
            title_string = 'Volumetric soil moisture (fraction) valid at ' \
              + '\n' + valid_when + ' ' \
              + str(int(round(lvl_bounds[lvl_idx-1]*100))) \
              + '-' + str(int(round(lvl_bounds[lvl_idx]*100))) + ' cm' \
              + ' from LAPS'
        p.title(title_string)
    return 

def plot_soilt(file_name, cntr_lvl=None):
    file, vars = peek(file_name, show_vars=False)
    lon = vars['lon'].get_value()
    lat = vars['lat'].get_value()
    soilt_var = vars['soil_temp']
    soilt = soilt_var.get_value()[0]
    mask = n.zeros(soilt.shape, dtype=int)
    for lvl_idx in range(len(soilt)):
        mask[lvl_idx] = vars['sea_lnd_ice_mask'].get_value()
    soilt_m = ma.masked_where(mask == 0 , soilt)
    lvl_bounds = vars['soil_lvl'].get_value()
    valid_date = str(vars['valid_date'].get_value()[0])
    valid_time = zfill(str(vars['valid_time'].get_value()[0]), 4)
    valid_when = valid_date[6:] + ' ' \
      + cardinal_2_month(int(valid_date[4:6])) + ' ' \
      + valid_date[0:4] \
      + ' ' + valid_time[:2] + ':' \
      + valid_time[2:] + ' UTC'
    m = set_default_basemap(lon,lat)
    # must plot using 2d lon and lat
    LON, LAT = p.meshgrid(lon,lat)
    #for lvl_idx in range(1):
    for lvl_idx in range(len(soilt)):
        p.figure()
        if cntr_lvl is not None:
            m.contourf(LON,LAT,soilt_m[lvl_idx], cntr_lvl)
        else:
            m.contourf(LON,LAT,soilt_m[lvl_idx])
        m.drawcoastlines()
        m.drawmeridians(n.array(n.arange(lon.min(), lon.max() + a_small_number, 15.)), labels=[1,0,0,1])
        m.drawparallels(n.array(n.arange(lat.min(), lat.max() + a_small_number, 15.)), labels=[1,0,0,1])
        p.colorbar(orientation='horizontal', shrink=0.7, fraction=0.02, pad=0.07, aspect=50)
        if lvl_idx == 0:
            title_string = 'Soil temperature (K) valid at ' \
              + '\n' + valid_when + ' ' \
              + '0' \
              + '-' + str(int(round(lvl_bounds[lvl_idx]*100))) + ' cm' \
              + ' from LAPS'
        else:
            title_string = 'Soil temperature (K) valid at ' \
              + '\n' + valid_when + ' ' \
              + str(int(round(lvl_bounds[lvl_idx-1]*100))) \
              + '-' + str(int(round(lvl_bounds[lvl_idx]*100))) + ' cm' \
              + ' from LAPS'
        p.title(title_string)
    return 

def plot_temp(file_name, cntr_lvl=None, save_frames=False):
    file, vars = peek(file_name, show_vars=False)
    lon = vars['lon'].get_value()
    lat = vars['lat'].get_value()
    air_temp_var = vars['air_temp']
    air_temp = air_temp_var.get_value()[0]
    lvl = vars['lvl'].get_value()
    valid_date = str(vars['valid_date'].get_value()[0])
    valid_time = zfill(str(vars['valid_time'].get_value()[0]), 4)
    valid_when = valid_date[6:] + ' ' \
      + cardinal_2_month(int(valid_date[4:6])) + ' ' \
      + valid_date[0:4] \
      + ' ' + valid_time[:2] + ':' \
      + valid_time[2:] + ' UTC'
    if save_frames:
        frame_number = 0
    m = set_default_basemap(lon,lat)
    # must plot using 2d lon and lat
    LON, LAT = p.meshgrid(lon,lat)
    #for lvl_idx in range(1):
    for lvl_idx in range(len(lvl)):
        p.figure()
        if cntr_lvl is not None:
            #p.contourf(lon,lat,air_temp_m[lvl_idx], cntr_lvl)
            #p.contourf(lon,lat,air_temp[lvl_idx], cntr_lvl)
            m.contourf(LON,LAT,air_temp[lvl_idx], cntr_lvl)
        else:
            #p.contourf(lon,lat,air_temp_m[lvl_idx])
            #p.contourf(lon,lat,air_temp[lvl_idx])
            m.contourf(LON,LAT,air_temp[lvl_idx])
        m.drawcoastlines()
        m.drawmeridians(n.array(n.arange(lon.min(), lon.max() + a_small_number, 15.)), labels=[1,0,0,1])
        m.drawparallels(n.array(n.arange(lat.min(), lat.max() + a_small_number, 15.)), labels=[1,0,0,1])
        p.colorbar(orientation='horizontal', shrink=0.7, fraction=0.02, pad=0.07, aspect=50)
        title_string = 'Air temperature (K) valid at ' \
          + '\n' + valid_when + ' ' \
          + str(int(round(lvl[lvl_idx]))) + ' mb' \
          + ' from LAPS'
        p.title(title_string)
        if save_frames:
            p.savefig('frames/frame_' + zfill(str(frame_number),3) +'_air_temp_' + str(int(lvl[lvl_idx])) + '_mb.png')
            frame_number += 1
    return 

# plot_wind
# plot_mr
def plot_mr(file_name, cntr_lvl=None, save_frames=False):
    file, vars = peek(file_name, show_vars=False)
    lon = vars['lon'].get_value()
    lat = vars['lat'].get_value()
    mr_var = vars['mix_rto']
    # to change between kg/kg to g/kg
    mr = mr_var.get_value()[0]*1000.
    lvl = vars['lvl'].get_value()
    valid_date = str(vars['valid_date'].get_value()[0])
    valid_time = zfill(str(vars['valid_time'].get_value()[0]), 4)
    valid_when = valid_date[6:] + ' ' \
      + cardinal_2_month(int(valid_date[4:6])) + ' ' \
      + valid_date[0:4] \
      + ' ' + valid_time[:2] + ':' \
      + valid_time[2:] + ' UTC'
    if save_frames:
        frame_number = 0
    m = set_default_basemap(lon,lat)
    # must plot using 2d lon and lat
    LON, LAT = p.meshgrid(lon,lat)
    #for lvl_idx in range(2):
    for lvl_idx in range(len(lvl)):
        p.figure()
        if cntr_lvl is not None:
            m.contourf(LON,LAT,mr[lvl_idx], cntr_lvl)
        else:
            m.contourf(LON,LAT,mr[lvl_idx])
        m.drawcoastlines()
        m.drawmeridians(n.array(n.arange(lon.min(), lon.max() + a_small_number, 15.)), labels=[1,0,0,1])
        m.drawparallels(n.array(n.arange(lat.min(), lat.max() + a_small_number, 15.)), labels=[1,0,0,1])
        p.colorbar(orientation='horizontal', shrink=0.7, fraction=0.02, pad=0.07, aspect=50)
        title_string = 'Mixing ration (g/kg)' \
          + '\n' + valid_when + ' ' \
          + str(int(round(lvl[lvl_idx]))) + ' mb' \
          + ' from LAPS'
        p.title(title_string)
        if save_frames:
            p.savefig('frames/frame_' + zfill(str(frame_number),3) +'_mr_' + str(int(lvl[lvl_idx])) + '_mb.png')
            frame_number += 1
    return 

# plot_ght
def plot_ght(file_name, cntr_lvl=None, save_frames=False):
    file, vars = peek(file_name, show_vars=False)
    lon = vars['lon'].get_value()
    lat = vars['lat'].get_value()
    ght_var = vars['geop_ht']
    ght = ght_var.get_value()[0]
    lvl = vars['lvl'].get_value()
    valid_date = str(vars['valid_date'].get_value()[0])
    valid_time = zfill(str(vars['valid_time'].get_value()[0]), 4)
    valid_when = valid_date[6:] + ' ' \
      + cardinal_2_month(int(valid_date[4:6])) + ' ' \
      + valid_date[0:4] \
      + ' ' + valid_time[:2] + ':' \
      + valid_time[2:] + ' UTC'
    if save_frames:
        frame_number = 0
    m = set_default_basemap(lon,lat)
    # must plot using 2d lon and lat
    LON, LAT = p.meshgrid(lon,lat)
    #for lvl_idx in range(2):
    for lvl_idx in range(len(lvl)):
        p.figure()
        if cntr_lvl is not None:
            m.contourf(LON,LAT,ght[lvl_idx], cntr_lvl)
        else:
            m.contourf(LON,LAT,ght[lvl_idx])
        m.drawcoastlines()
        m.drawmeridians(n.array(n.arange(lon.min(), lon.max() + a_small_number, 15.)), labels=[1,0,0,1])
        m.drawparallels(n.array(n.arange(lat.min(), lat.max() + a_small_number, 15.)), labels=[1,0,0,1])
        p.colorbar(orientation='horizontal', shrink=0.7, fraction=0.02, pad=0.07, aspect=70)
        title_string = 'Geopotential height (m)' \
          + '\n' + valid_when + ' ' \
          + str(int(round(lvl[lvl_idx]))) + ' mb' \
          + ' from LAPS'
        p.title(title_string)
        if save_frames:
            p.savefig('frames/frame_' + zfill(str(frame_number),3) +'_ght_' + str(int(lvl[lvl_idx])) + '_mb.png')
            frame_number += 1
    return 

def plot_sfc_temp(file_name, cntr_lvl=None, save_frames=False):
    file, vars = peek(file_name, show_vars=False)
    lon = vars['lon'].get_value()
    lat = vars['lat'].get_value()
    sfc_temp_var = vars['sfc_temp']
    sfc_temp = sfc_temp_var.get_value()[0]
    valid_date = str(vars['valid_date'].get_value()[0])
    valid_time = zfill(str(vars['valid_time'].get_value()[0]), 4)
    valid_when = valid_date[6:] + ' ' \
      + cardinal_2_month(int(valid_date[4:6])) + ' ' \
      + valid_date[0:4] \
      + ' ' + valid_time[:2] + ':' \
      + valid_time[2:] + ' UTC'
    m = set_default_basemap(lon,lat)
    # must plot using 2d lon and lat
    LON, LAT = p.meshgrid(lon,lat)
    p.figure()
    if cntr_lvl is not None:
        m.contourf(LON,LAT,sfc_temp, cntr_lvl)
    else:
        m.contourf(LON,LAT,sfc_temp)
    m.drawcoastlines()
    m.drawmeridians(n.array(n.arange(lon.min(), lon.max() + a_small_number, 15.)), labels=[1,0,0,1])
    m.drawparallels(n.array(n.arange(lat.min(), lat.max() + a_small_number, 15.)), labels=[1,0,0,1])
    p.colorbar(orientation='horizontal', shrink=0.7, fraction=0.02, pad=0.07, aspect=70)
    title_string = 'Surface temperature (K) valid at' \
      + '\n' + valid_when + ' ' \
      + ' from LAPS'
    p.title(title_string)
    if save_frames:
        p.savefig('frames/frame_' + zfill(str(frame_number),3) +'_sfc_temp_' + str(int(lvl[lvl_idx])) + '.png')
    return 

def plot_sfc_pres(file_name, cntr_lvl=None, save_frames=False):
    file, vars = peek(file_name, show_vars=False)
    lon = vars['lon'].get_value()
    lat = vars['lat'].get_value()
    sfc_pres_var = vars['sfc_pres']
    sfc_pres = sfc_pres_var.get_value()[0]/100.
    valid_date = str(vars['valid_date'].get_value()[0])
    valid_time = zfill(str(vars['valid_time'].get_value()[0]), 4)
    valid_when = valid_date[6:] + ' ' \
      + cardinal_2_month(int(valid_date[4:6])) + ' ' \
      + valid_date[0:4] \
      + ' ' + valid_time[:2] + ':' \
      + valid_time[2:] + ' UTC'
    m = set_default_basemap(lon,lat)
    # must plot using 2d lon and lat
    LON, LAT = p.meshgrid(lon,lat)
    p.figure()
    if cntr_lvl is not None:
        m.contourf(LON,LAT,sfc_pres, cntr_lvl)
    else:
        m.contourf(LON,LAT,sfc_pres)
    m.drawcoastlines()
    m.drawmeridians(n.array(n.arange(lon.min(), lon.max() + a_small_number, 15.)), labels=[1,0,0,1])
    m.drawparallels(n.array(n.arange(lat.min(), lat.max() + a_small_number, 15.)), labels=[1,0,0,1])
    p.colorbar(orientation='horizontal', shrink=0.7, fraction=0.02, pad=0.07, aspect=70)
    title_string = 'Surface pressure (K) valid at' \
      + '\n' + valid_when + ' ' \
      + ' from LAPS'
    p.title(title_string)
    if save_frames:
        p.savefig('frames/frame_' + zfill(str(0),3) +'_sfc_pres.png')
    return 

# plot_sfc_wind
# plot_mslp
# plot_land
# plot_terrain
def plot_skin_temp(file_name, cntr_lvl=None, save_frames=False):
    file, vars = peek(file_name, show_vars=False)
    lon = vars['lon'].get_value()
    lat = vars['lat'].get_value()
    skin_temp = vars['skin_temp']
    skin_temp = skin_temp_var.get_value()[0]
    valid_date = str(vars['valid_date'].get_value()[0])
    valid_time = zfill(str(vars['valid_time'].get_value()[0]), 4)
    valid_when = valid_date[6:] + ' ' \
      + cardinal_2_month(int(valid_date[4:6])) + ' ' \
      + valid_date[0:4] \
      + ' ' + valid_time[:2] + ':' \
      + valid_time[2:] + ' UTC'
    m = set_default_basemap(lon,lat)
    # must plot using 2d lon and lat
    LON, LAT = p.meshgrid(lon,lat)
    p.figure()
    if cntr_lvl is not None:
        m.contourf(LON,LAT,skin_temp, cntr_lvl)
    else:
        m.contourf(LON,LAT,skin_temp)
    m.drawcoastlines()
    m.drawmeridians(n.array(n.arange(lon.min(), lon.max() + a_small_number, 15.)), labels=[1,0,0,1])
    m.drawparallels(n.array(n.arange(lat.min(), lat.max() + a_small_number, 15.)), labels=[1,0,0,1])
    p.colorbar(orientation='horizontal', shrink=0.7, fraction=0.02, pad=0.07, aspect=70)
    title_string = 'Surface pressure (K) valid at' \
      + '\n' + valid_when + ' ' \
      + ' from LAPS'
    p.title(title_string)
    if save_frames:
        p.savefig('frames/frame_' + zfill(str(frame_number),3) +'_skin_temp_' + str(int(lvl[lvl_idx])) + '.png')
    return 

# plot_skin_temp


