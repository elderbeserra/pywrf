import sys
sys.path.append('/Users/val/data/pylib')

import viz.utils as vu
from string import zfill

data_file = \
  '/Volumes/DATA/wrf/first_1way/coarse/wrfout_d01_2003-01-17_03-00-00.nc'
f, fv = vu.peek(data_file, show_vars=False)

u = fv['U'].get_value()
lon = fv['XLONG_U'].get_value()
lat = fv['XLAT_U'].get_value()

lvl_idx = 10
times = fv['Times'].get_value()

for time_idx in range(46):
    time_string = vu.set_time_string('WRF', times[time_idx])
    lvl_string = str(lvl_idx) + 'th lvl'
    title_string = \
      vu.set_title_string('u wind', 'm/s', time_string, lvl_string, 'WRF:')
    file_name = 'frames/frame_' + zfill(str(time_idx),3)
    print file_name
    vu.plot_slab(lon[time_idx], lat[time_idx], u[time_idx][lvl_idx], 
      file_name=file_name, title_string=title_string)
