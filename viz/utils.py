"""
repo of useful functions for visualizing stuff
"""

import PyNGL_numpy.Nio as nio
import pylab as p
import matplotlib.numerix.ma as ma
import numpy as n
from string import zfill
from matplotlib.toolkits.basemap import Basemap

a_small_number = 1e-8

def get_long_name(var):
    try: 
        long_name = var.long_name
    except AttributeError:
        try:
            long_name = var.description
        except AttributeError:
            long_name = 'N/A'
    if long_name == '':
        long_name = 'N/A'
    return long_name

def peek(file_name, return_pointers=False, show_vars=True):
    #print file_name
    file = nio.open_file(file_name)
    vars = file.variables
    if show_vars:
        keys = vars.keys()
        keys.sort()
        page_idx = 0
        max = 0
        lines_in_a_page = 50
        if len(keys) > lines_in_a_page:
            while max < len(keys):
                min = page_idx * lines_in_a_page
                max = (page_idx + 1) * lines_in_a_page
                if max > len(keys):
                    max = len(keys)
                print min, max, type(min), type(max)
                for key in keys[min:max]:
                    long_name = get_long_name(vars[key])
                    if len(key) < 7:
                        print '\t', key, '\t\t->\t', long_name
                    else:
                        print '\t', key, '\t->\t', long_name
                page_idx += 1
                raw_input('press enter to continue...')
        else:
            for key in keys:
                long_name = get_long_name(vars[key])
                if len(key) < 7:
                    print '\t', key, '\t\t->\t', long_name
                else:
                    print '\t', key, '\t->\t', long_name
    if return_pointers:
        return file, vars
    else:
        return

def cardinal_2_month(cardinal):
    if cardinal == 1:
        return 'Jan'
    elif cardinal == 2:
        return 'Feb'
    elif cardinal == 3:
        return 'Mar'
    elif cardinal == 4:
        return 'Apr'
    elif cardinal == 5:
        return 'May'
    elif cardinal == 6:
        return 'Jun'
    elif cardinal == 7:
        return 'Jul'
    elif cardinal == 8:
        return 'Aug'
    elif cardinal == 9:
        return 'Sep'
    elif cardinal == 10:
        return 'Oct'
    elif cardinal == 11:
        return 'Nov'
    elif cardinal == 12:
        return 'Dec'
    else:
        return 'gibberish'

def set_default_basemap(lon, lat, frame_width=5.):
    test = lon < 0.
    if True in test:
        # matplotlib expects 0-360 while WRF for example uses -180-180
        delta = n.ones(lon.shape)
        delta *= 360
        delta = ma.masked_where(lon > 0., delta)
        lon += delta.filled(fill_value=0)
    llcrnrlon=lon.min() - frame_width
    urcrnrlon=lon.max() + frame_width
    llcrnrlat=lat.min() - frame_width
    urcrnrlat=lat.max() + frame_width
    lon_0 = llcrnrlon + (urcrnrlon - llcrnrlon) / 2.
    lat_0 = llcrnrlat + (urcrnrlat - llcrnrlat) / 2.
        
    map = Basemap(
      llcrnrlon=llcrnrlon,
      llcrnrlat=llcrnrlat,
      urcrnrlon=urcrnrlon,
      urcrnrlat=urcrnrlat,
      resolution='l',
      projection='cyl',
      lon_0=lon_0,
      lat_0=lat_0
      )
    return map

def set_lvl_string(lvl, lvl_type):
    if lvl_type == 'sigma':
        pass
    if lvl_type == 'press':
        pass
    if lvl_type == 'soil':
        pass
    return lvl_string

def set_title_string(
  long_name, 
  units,
  time_string, 
  lvl_string, 
  prefix='', 
  postfix=''):
    title_string = \
      prefix + ' '\
      + lvl_string + ' ' \
      + long_name + ' (' + units + ')' \
      + ' valid at\n' \
      + time_string \
      + postfix
    return title_string

def set_time_string(model, time):
    if model == 'WRF':
        time_string = time.tostring()
        year = time_string[:4]
        month = cardinal_2_month(int(time_string[5:7]))
        day = time_string[8:10]
        hour = time_string[11:13]
        minute = time_string[14:16]
        second = time_string[17:]
        time_string = day + ' ' \
          + month + ' ' \
          + year + ' ' \
          + hour + ':' + minute + ' UTC'
        return time_string
    elif model == 'LAPS':
        valid_date = str(time[0])
        valid_time = zfill(str(time[1]),4)
        time_string = valid_date[6:] + ' ' \
          + cardinal_2_month(int(valid_date[4:6])) + ' ' \
          + valid_date[0:4] \
          + ' ' + valid_time[:2] + ':' \
          + valid_time[2:] + ' UTC'
        return time_string
    elif model == 'manual':
        # expects list of the form
        # [year, month, day, hour, minute]
        time_string = (str(time[2]) + ' ')
        time_string += (cardinal_2_month(time[1]) + ' ')
        time_string += (str(time[0]) + ' ')
        time_string += (zfill(str(time[3]),2) + ':')
        time_string += (zfill(str(time[4]),2) + ' ')
        time_string += 'UTC'
        return time_string
    else:
        print 'set_time_string has done nothing'
        return

def plot_grid(lon,lat,
  title_string = 'N/A',
  meridians_delta = 15,
  parallels_delta = 15,
  same_figure = False,
  figsize = None,
  file_name = None,
  dpi = 80,
  skip = 5,
  return_map = False,
  marker = '+'
  ):
    """Function to plot grids similar to those generated by WPS
    
    Modified 27/01/08: Minor mod to plot the boundaries of the domain
    Comments: There is a chunk of code that gets repeated in a lot of places... the 
    stuff about plotting meridians and parallels. Maybe we could put this in a seperate 
    funtion?
    
    """
    map = set_default_basemap(lon,lat)
    # let's make sure the lat and lon arrays are 2D
    if len(lon.shape) < 2:
        lon, lat = p.meshgrid(lon,lat)
    if not same_figure:
        if not figsize:
            p.figure()
        else:
            p.figure(figsize=figsize)
    map.plot(lon[::skip,::skip], lat[::skip,::skip], marker=marker, linestyle='None')

    corners_lon = n.array([lon[0,0], lon[0,-1], lon[-1,-1], lon[-1,0]])
    corners_lat = n.array([lat[0,0], lat[0,-1], lat[-1,-1], lat[-1,0]])
    
    map.plot(corners_lon,corners_lat, 'ro')

    # Here it is =)
    map.plot(lon[0,:],lat[0,:],'k-')
    map.plot(lon[:,-1],lat[:,-1],'k-')
    map.plot(lon[-1,:],lat[-1,:],'k-')
    map.plot(lon[:,0],lat[:,0],'k-')

#    Thom: I've taken out the plot canberra option for generality
#    canberra_lon = [149 + 8./60]
#    canberra_lat = [-35 - 17./60]
#    map.plot(canberra_lon,canberra_lat, 'gs')

    if not same_figure:
        map.drawcoastlines()
	# These blocks seem to get repeated...
        meridians = n.arange(int(round(lon.min(),0)), 
          lon.max() + a_small_number, meridians_delta)
        meridians = n.array(meridians)
        map.drawmeridians(meridians, labels=[1,0,0,1])
    
        parallels = n.arange(int(round(lat.min(),0)), 
          lat.max() + a_small_number, parallels_delta)
        parallels = n.array(parallels)
        map.drawparallels(parallels, labels=[1,0,0,1])

        p.title(title_string)
    if file_name:
        p.savefig(file_name,dpi=dpi)
    if return_map:
        return map
    else:
        return

def plot_slab(lon, lat, slab,
  map = None,
  figsize = None,
  cntr_lvl = None,
  title_string = 'N/A',
  meridians_delta = 15,
  parallels_delta = 15,
  file_name = None,
  dpi = 80,
  wind_vector =None,
  quiv_skip = 5,
  frame_width = 5,
  significant_digits = 0,
  return_map = False
  ):
    # let's make sure the lat and lon arrays are 2D
    if len(lon.shape) < 2:
        lon, lat = p.meshgrid(lon,lat)
    if map == None:
        map = set_default_basemap(lon,lat,frame_width)
    if not figsize:
        p.figure()
    else:
        p.figure(figsize=figsize)
    if cntr_lvl != None:
        map.contourf(lon, lat, slab, cntr_lvl)
    else:
        map.contourf(lon, lat, slab)
    if wind_vector != None:
        map.quiver(lon[::quiv_skip,::quiv_skip], 
          lat[::quiv_skip,::quiv_skip],
          wind_vector[0][::quiv_skip,::quiv_skip],
          wind_vector[1][::quiv_skip,::quiv_skip])

    map.drawcoastlines()

    # todo
    # the +5 is a hack in case one uses the default map it should be made an
    # an explicit parameter that is passed around...
    meridians = n.arange(int(round(lon.min(),0)), 
      #int(round(lon.max(),0)) + a_small_number, meridians_delta)
      int(round(lon.max(),0)) + 5 + a_small_number, meridians_delta)
    meridians = n.array(meridians)
    map.drawmeridians(meridians, labels=[1,0,0,1])

    parallels = n.arange(int(round(lat.min(),0)), 
      #int(round(lat.max(),0)) + a_small_number, parallels_delta)
      int(round(lat.max(),0)) + 5 + a_small_number, parallels_delta)
    parallels = n.array(parallels)
    map.drawparallels(parallels, labels=[1,0,0,1])

    #plot canberra
    if 1:
        map.plot([149.133],[-35.283],'bo', ms=3)

    format = '%.'+ str(significant_digits) + 'f'
    p.colorbar(orientation='horizontal', shrink=0.7, 
      fraction=0.02, pad=0.07, aspect=70, 
      format = format)
    p.title(title_string)
    if file_name:
        p.savefig(file_name,dpi=dpi)
        p.close()
    if return_map:
        return map
    else:
        return

def flip_yaxis(slab):
    shape = slab.shape
    # assuming (lat, lon) ordering
    lat_idx_max = shape[0] - 1
    flipped_slab = n.zeros(shape)
    for lat_idx in range(lat_idx_max, -1, -1):
        flipped_idx = abs(lat_idx - lat_idx_max)
        flipped_slab[flipped_idx] = slab[lat_idx]
    return flipped_slab

 
def set_cntr_lvl(data, intervals=12, only_positive=False, significant_digits=0):
    import math as m
    dummy = data * 10**significant_digits
    min = dummy.min()
    max = dummy.max()
    floor = m.floor(min)
    ceil = m.ceil(max)
    if only_positive:
        floor = 0
    extent = float(ceil - floor)
    delta = extent / intervals
    #cntr_lvl = n.arange(floor, ceil + a_small_number, delta)
    cntr_lvl = n.arange(floor - a_small_number, ceil + a_small_number, delta)
    return cntr_lvl / 10**significant_digits

def lin_interp(xa, xb, ya, yb, x):
    """linear interpolation in two dimensions
    
    USAGE
    xa = original x-grid
    ya = original y-grid
    xb = new x-grid
    yb = new y-grid
    x  = data in xa, ya

    """
    slope = (yb - ya) / (xb - xa)
    y = ya + slope * (x - xa)
    return y
    



