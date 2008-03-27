"""
repo of useful functions for visualizing stuff
"""

import os
from socket import gethostname

# VB The following import is host specific
hostname = gethostname()
if hostname == 'hn3.its.monash.edu.au':
    import PyNGL.Nio as nio
elif hostname == 'linux450':
    # VB Sorry Thom if this is not correct ;)
    import PyNGL.Nio as nio
elif hostname == 'val.maths.monash.edu.au':
    import Nio as nio
    #import PyNGL_numpy.Nio as nio
elif hostname == 'valerio-bisignanesis-computer.local':
    import Nio as nio
else:
    print 'Warning: since I do not know this hostname, I am not sure of ' \
      + ' the appropriate syntax to import pynio... I will try\n' \
      + ' import PyNGL.Nio as nio'
    import PyNGL.Nio as nio

import time
import pylab as p
import matplotlib.numerix.ma as ma
import numpy as n
from string import zfill
from matplotlib.toolkits.basemap import Basemap
import gc

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
      + long_name + ' (' + units + ')\n' \
      + ' valid at ' \
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
  colorbar = False,
  contour_labels = True,
  monochrome = False,
  return_map = False
  ):
    from matplotlib.cm import gist_ncar as cmap
    # let's make sure the lat and lon arrays are 2D
    if len(lon.shape) < 2:
        lon, lat = p.meshgrid(lon,lat)
    if map is None:
        map = set_default_basemap(lon,lat,frame_width)
    if not figsize:
        fig = p.figure()
    else:
        fig = p.figure(figsize=figsize)
    if cntr_lvl is not None:
        if not monochrome:
            csetf = map.contourf(lon, lat, slab, 
              cntr_lvl, 
              cmap=cmap)
        cset = map.contour(lon, lat, slab, cntr_lvl, colors='lightslategray')
    else:
        if not monochrome:
            csetf = map.contourf(lon, lat, slab, cmap=cmap)
        cset = map.contour(lon, lat, slab, colors='lightslategray')
    if wind_vector is not None:
        quiv = map.quiver(lon[::quiv_skip,::quiv_skip], 
          lat[::quiv_skip,::quiv_skip],
          wind_vector[0][::quiv_skip,::quiv_skip],
          wind_vector[1][::quiv_skip,::quiv_skip])
    # plot grid outline
    map.plot(lon[0,:],lat[0,:],color='lightslategray')
    map.plot(lon[-1,:],lat[-1,:],color='lightslategray')
    map.plot(lon[:,0],lat[:,0],color='lightslategray')
    map.plot(lon[:,-1],lat[:,-1],color='lightslategray')
    if monochrome:
        map.fillcontinents(color='0.95')

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

    if contour_labels:
        p.clabel(cset, cset.levels[::2], 
          colors='k', fontsize=8, fmt='%i')

    if colorbar:
        format = '%.'+ str(significant_digits) + 'f'
        p.colorbar(csetf, orientation='horizontal', shrink=0.7, 
          fraction=0.02, pad=0.07, aspect=70, 
          format = format)

    p.title(title_string)
    if file_name:
        p.savefig(file_name,dpi=dpi)
        p.close(fig)
        del fig
    if not monochrome:
        del csetf
    del cset
    if wind_vector is not None:
        del quiv
    gc.collect()
    if return_map:
        return map
    else:
        del map
        gc.collect()
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

 
def set_cntr_lvl(data, 
      intervals=12, 
      only_positive=False, 
      significant_digits=0,
      forced_max=None,
      forced_min=None):

    import math as m
    dummy = data * 10**significant_digits
    min = dummy.min()
    max = dummy.max()
    floor = m.floor(min)
    ceil = m.ceil(max)
    if forced_max:
        ceil = forced_max * 10**significant_digits
    if forced_min:
        floor = forced_min * 10**significant_digits
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
    
def plot_slice(
  ordinate,
  data,
  abscissa = None,
  abscissa_label = None,
  land_mask = None,
  cntr_lvl = None,
  wind_vector = None,
  log_scale = False,
  pressure = True,
  ordinate_quiv_skip = 1,
  abscissa_quiv_skip = 3,
  xmin = None,
  xmax = None,
  ymin = None,
  ymax = None,
  significant_digits = 0,
  title_string = 'N/A',
  file_name = None,
  dpi = 100
  ):
    '''
    plot a slice
    Usage:
    >>> plot_vertical_slice(ordinate, data, abscissa, wind_vector)
    where
      ordinate is either pressure or height. By default pressure is assumed
      data is a vertical slice
      abscissa is either 1D or 2D
    '''
    if log_scale:
        ordinate = n.log10(ordinate)
    # if the abscissa is not supplied than simply use the record numbers
    if abscissa is None:
        x = n.arange(1,data.shape[1]+1)
        abscissa = n.zeros(data.shape)
        for y_idx in range(data.shape[0]):
            abscissa[y_idx] = x
        if cntr_lvl is not None:
            cset = p.contourf(abscissa, ordinate, data, cntr_lvl)
        else:
            cset = p.contourf(abscissa, ordinate, data)
        p.xlabel('record number')
    else:
        # let's handle 1D abscissa arrays
        if len(abscissa.shape) == 1:
            dummy = n.zeros(data.shape)
            for lvl_idx in range(data.shape[0]):
                dummy[lvl_idx] = abscissa
            abscissa = dummy
            del dummy
        if cntr_lvl is not None:
            cset = p.contourf(abscissa, ordinate, data, cntr_lvl)
        else:
            cset = p.contourf(abscissa, ordinate, data)
        if abscissa_label:
            p.xlabel(abscissa_label)
        else:
            p.xlabel('N/A')
    if wind_vector:
        p.quiver(abscissa[::ordinate_quiv_skip,::abscissa_quiv_skip], 
          ordinate[::ordinate_quiv_skip,::abscissa_quiv_skip], 
          wind_vector[0][::ordinate_quiv_skip,::abscissa_quiv_skip], 
          wind_vector[1][::ordinate_quiv_skip,::abscissa_quiv_skip])
    if land_mask is not None:
        land = ma.masked_where(land_mask == 2,ordinate[0])
        p.plot(abscissa[0], land, color=(0.59,0.29,0.0), linewidth=2.)
        # if you also want to plot the ocean uncomment the following lines
        #if log_scale:
        #    ocean = ma.masked_where(land_mask == 1, n.log10(1013.25))
        #else:
        #    ocean = ma.masked_where(land_mask == 1, 1013.25)
        #p.plot(abscissa[0], ocean, color=(0.0,0.0,1.0), linewidth=2.)
    if pressure:
        # I am assuming pressure will be expressed in hPa
        yticks_location = n.arange(1000.,99.,-100.)
        if log_scale:
            p.yticks(n.log10(yticks_location), 
              [str(int(e)) for e in yticks_location])
            p.ylabel('log10 of pressure')
            for e in n.log10(yticks_location):
                p.axhline(e,linestyle='--',color=(0.7,0.7,0.7))
            # the -5. is there to create a bit of a top margin
            if ordinate.max() > n.log10(1013.25):
                cheat_ymin = ordinate.max()
            else:
                cheat_ymin = n.log10(1013.25 + 5.)
            p.ylim(ymin=cheat_ymin,
              ymax=n.log10(10**ordinate.min() - 5.))
        else:
            p.yticks( yticks_location, 
              [str(int(e)) for e in yticks_location])
            p.ylabel('pressure (hPa)')
            for e in yticks_location:
                p.axhline(e,linestyle='--',color=(0.7,0.7,0.7))
            # the -25. is there to create a bit of a top margin
            if ordinate.max() > 1013.25:
                cheat_ymin = ordinate.max()
            else:
                cheat_ymin = 1013.25 + 10.
            p.ylim(ymin=cheat_ymin, ymax=ordinate.min() - 25.)
            # p.ylim(ymin=ordinate.max(), ymax=ordinate.min() - 25.)
            # if any of the axes boundaries have been given explicitly, we'll 
            # them
            if log_scale:
                if xmin is not None:
                    p.xlim(xmin=n.log10(xmin))
                if xmax is not None:
                    p.xlim(xmax=n.log10(xmax))
                if ymin is not None:
                    p.ylim(ymin=n.log10(ymin))
                if ymax is not None:
                    p.ylim(ymax=n.log10(ymax))
            else:
                if xmin is not None:
                    p.xlim(xmin=xmin)
                if xmax is not None:
                    p.xlim(xmax=xmax)
                if ymin is not None:
                    p.ylim(ymin=ymin)
                if ymax is not None:
                    p.ylim(ymax=ymax)

    else:
        print 'I assume you are trying to plot in z coordinate: ' \
          + 'sorry not implemented yet'
    format = '%.'+ str(significant_digits) + 'f'
    p.colorbar(cset, orientation='horizontal', shrink=0.7, 
      #fraction=0.02, pad=0.095, aspect=70, 
      fraction=0.02, pad=0.1, aspect=70, 
      format = format)
    p.title(title_string)
    if file_name:
        p.savefig(file_name,dpi=dpi)
        p.close()

def write_to_log_file(log_file, message):
    '''
    This functions opens, writes to it,
    and closes the log file so that tools like tail -f
    will work as intended.
    It also prepends a time stamp to each entry.
    It is assumed that the output_dir and log_file are defined in the 
    namespace from which the function is called.
    '''
    log_file = open(log_file, 'a')
    log_file.write(time.ctime(time.time()) + ' -> ' + message + '\n')
    log_file.close()

def generate_output_file_name(output_dir, prefix, timetuple):
    """Returns the output file name built by joining the output directory
    path with the supplied prefix and the timetuple which is used to
    construct the suffix/timestamp
    """
    output_file_name = prefix
    output_file_name += str(timetuple[0])
    output_file_name += str(timetuple[1]).zfill(2)
    output_file_name += str(timetuple[2]).zfill(2)
    output_file_name += str(timetuple[3]).zfill(2)
    output_file_name += str(timetuple[4]).zfill(2)
    output_file_name += str(timetuple[5]).zfill(2)
    output_file_name = os.path.join(output_dir, output_file_name)
    return output_file_name
    
