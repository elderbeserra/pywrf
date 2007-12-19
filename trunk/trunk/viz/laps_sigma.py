from viz.utils import peek
from viz.utils import cardinal_2_month
import pylab as p
import matplotlib.numerix.ma as ma
import numpy as n
from string import zfill
from matplotlib.toolkits.basemap import Basemap

a_small_number = 1e-8

def set_default_basemap(lon, lat):
    map = Basemap(
      llcrnrlon=lon.min() - 0.5,
      llcrnrlat=lat.min() - 0.5,
      urcrnrlon=lon.max() + 0.5,
      urcrnrlat=lat.max() + 0.5,
      resolution='l',
      projection='cyl',
      lon_0=lon.min() + (lon.max()-lon.min())/2.,
      lat_0=lat.min() + (lat.max()-lat.min())/2.
      )
    return map

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
    for lvl_idx in range(5):
    #for lvl_idx in range(len(lvl)):
        p.figure()
        if cntr_lvl != None:
            m.contourf(LON,LAT,air_temp[lvl_idx], cntr_lvl)
        else:
            m.contourf(LON,LAT,air_temp[lvl_idx])
        m.drawcoastlines()
        m.drawmeridians(n.array(n.arange(lon.min(), lon.max() + a_small_number, 15.)), labels=[1,0,0,1])
        m.drawparallels(n.array(n.arange(lat.min(), lat.max() + a_small_number, 15.)), labels=[1,0,0,1])
        p.colorbar(orientation='horizontal', shrink=0.7, fraction=0.02, pad=0.07, aspect=50)
        title_string = 'Air temperature (K) valid at ' \
          + '\n' + valid_when + ' ' \
          + str(round(lvl[lvl_idx],4)) + ' sigma' \
          + ' from LAPS'
        p.title(title_string)
        if save_frames:
            p.savefig('frames/frame_' + zfill(str(frame_number),3) +'_air_temp_' + str(lvl[lvl_idx]) + '_sigma.png')
            frame_number += 1
    return 


