import os
import pylab as p
import matplotlib.numerix.ma as ma
import numpy as n
from string import zfill
from matplotlib.toolkits.basemap import Basemap

a_small_number = 1e-8

def colons_to_underscores(test=False):
    files = os.listdir('.')
    print 'I plan to do the following:'
    for file in files:
        command = 'mv ' + file + ' ' + file[:-6] + '_00_00'
        print command
    print 'Happy? [y/n]'
    answer = raw_input()
    if answer == 'y':
        for file in files:
            command = 'mv ' + file + ' ' + file[:-6] + '_00_00'
            os.system(command)
    else:
        print 'OK, maybe next time ;]' 
    return

def wrf_to_pressure(var, pressure, press_lvl, fill_value=1e35, positive_only=False):
    # these are the pressure levels we want as output
    k_max, j_max, i_max = pressure.shape
    output = n.zeros((len(press_lvl), j_max, i_max))
    for lvl_idx in range(len(press_lvl)):
        for j in range(j_max):
            for i in range(i_max):
                # make it ascending for searchsorted
                pressure_column = pressure[::-1,j,i]
                if press_lvl[lvl_idx] > pressure_column.max() \
                  or press_lvl[lvl_idx] < pressure_column.min():
                    output[lvl_idx,j,i] = fill_value
                else:
                    var_column = var[::-1,j,i]
                    press_idx = pressure_column.searchsorted(press_lvl[lvl_idx])
                    xa = pressure_column[press_idx]
                    xb = pressure_column[press_idx - 1]
                    ya = var_column[press_idx]
                    yb = var_column[press_idx - 1]
                    x = press_lvl[lvl_idx]
                    y = lin_interp(xa, xb, ya, yb, x)
                    if positive_only and y < 0.:
                        y = 0.
                    output[lvl_idx,j,i] = y
    return output

def wrf2latlon(projection, 
  lat_0, lon_0,
  lat_1,
  lat_2,
  grid_centre_lon,
  grid_centre_lat,
  nx,
  ny,
  delta,
  staggered = False,
  return_extra = False
  ):
    from pyproj import Proj
    proj = Proj(proj=projection, lat_0=lat_0, lon_0=lon_0, lat_1=lat_1, lat_2=lat_2)
    grid_centre_x, grid_centre_y = proj(grid_centre_lon, grid_centre_lat)
    grid_x_extent = nx * delta
    grid_y_extent = ny * delta
    min_x = grid_centre_x - grid_x_extent/2.
    min_y = grid_centre_y - grid_y_extent/2.
    max_x = min_x + grid_x_extent
    max_y = min_y + grid_y_extent
    x = n.arange(min_x, max_x + a_small_number, delta)
    y = n.arange(min_y, max_y + a_small_number, delta)  
    X, Y = p.meshgrid(x,y)
    lon, lat = proj(X, Y, inverse=True)

    if staggered:
        x_u = n.arange(min_x, max_x + delta + a_small_number, delta)
        x_u -= (delta /2.)
        X_u, Y_u = p.meshgrid(x_u,y)
        y_v = n.arange(min_y, max_y + delta + a_small_number, delta)  
        y_v -= (delta /2.)
        X_v, Y_v = p.meshgrid(x,y_v)
        lon_u, lat_u = proj(X_u, Y_u, inverse=True)
        lon_v, lat_v = proj(X_v, Y_v, inverse=True)

    llcrnrlon, llcrnrlat = proj(min_x, min_y, inverse=True)
    urcrnrlon, urcrnrlat = proj(max_x, max_y, inverse=True)
    map = Basemap(
      projection = projection,
      lon_0 = lon_0,
      lat_0 = lat_0,
      lat_1 = lat_1,
      lat_2 = lat_2,
      llcrnrlon = llcrnrlon,
      llcrnrlat = llcrnrlat,
      urcrnrlon = urcrnrlon,
      urcrnrlat = urcrnrlat
      )

    if return_extra:
        # it seems that the basemap automatically sets the origin the native
        # coordinate system at its llcrnr (or close to it...)
        offset = map(lon_0, lat_0)
        if staggered:
            X += offset[0]
            Y += offset[1]
            X_u += offset[0]
            Y_u += offset[1]
            X_v += offset[0]
            Y_v += offset[1]
            return lon, lat, lon_u, lat_u, lon_v, lat_v, \
              X, Y, X_u, Y_u, X_v, Y_v, map
        else:
            X += offset[0]
            Y += offset[1]
            return lon, lat, X, Y, map
    else: 
        if staggered:
            return lon, lat, lon_u, lat_u, lon_v, lat_v
        else:
            return lon, lat
  
def wrf_grid(
  # WPS -> map_proj
  projection,
  # WPS -> truelat1
  lat_1,
  # WPS -> truelat2
  lat_2,
  # WPS -> stand_lon
  lon_0,
  # WPS -> ref_lat
  grid_centre_lat,
  # WPS -> ref_lon
  grid_centre_lon,
  delta_x,
  delta_y,
  # WPS -> e_we
  nx,
  # WPS -> e_sn
  ny,
  show_mass_grid = False,
  show_stag_grids = False,
  ):
    if lon_0 != grid_centre_lon:
        print 'not implemented yet -> see the source'
        return

    width   = nx * delta_x
    height  = ny * delta_y
    frame_x = 10 * delta_x
    frame_y = 10 * delta_y
    m = Basemap(
      projection = projection,
      lat_0 = grid_centre_lat,
      # this could be a bad assumption... because lon_0 and grid_centre_lon
      # need not be aligned, but at the same time I need to give this to
      # basemap for the grid to be centred... I could probably fix it
      # assigning lon_0 and then imposing a grid shift in native coordinates
      # if ref_lon and lon_0 were not the same
      lon_0 = lon_0,
      lat_1 = lat_1,
      lat_2 = lat_2,
      width = width + 2*frame_x,
      height = height + 2*frame_y,
      resolution = 'l',
      area_thresh=1000.
      )
    grid_centre_x, grid_centre_y = m(grid_centre_lon, grid_centre_lat)
    min_x = grid_centre_x - width/2.
    min_y = grid_centre_y - height/2.
    max_x = min_x + width
    max_y = min_y + height
    x = n.arange(min_x, max_x + a_small_number, delta_x)
    y = n.arange(min_y, max_y + a_small_number, delta_y)  
    x = x[1:-1]
    y = y[1:-1]
    x_u = n.arange(min_x, max_x + delta_x + a_small_number, delta_x)
    x_u -= delta_x/2.
    x_u = x_u[1:-1]
    y_v = n.arange(min_y, max_y + delta_y + a_small_number, delta_y)
    y_v -= delta_y/2.
    y_v = y_v[1:-1]
    X, Y = p.meshgrid(x,y)
    lon, lat = m(X, Y, inverse=True)
    X_u, Y_u = p.meshgrid(x_u,y)
    lon_u, lat_u = m(X_u, Y_u, inverse=True)
    X_v, Y_v = p.meshgrid(x,y_v)
    lon_v, lat_v = m(X_v, Y_v, inverse=True)
    if show_mass_grid:
        m.plot(X, Y, 'b+')
        m.plot([grid_centre_x], [grid_centre_y], 'r+')
        if show_stag_grids:
            m.plot(X_u, Y_u, 'g+')
            m.plot(X_v, Y_v, 'r+')
        m.drawcoastlines()
	p.show()
    output = {
      'map' : m,
      'mass_stag': {
        'lon_2d' : lon,
        'lat_2d' : lat,
        'x'      : x,
        'y'      : y,
        'x_2d'   : X,
        'y_2d'   : Y,
        }, 
      'u_stag': {
        'lon_2d' : lon_u,
        'lat_2d' : lat_u,
        'x'      : x_u,
        'y'      : y,
        'x_2d'   : X_u,
        'y_2d'   : Y_u,
        }, 
      'v_stag': {
        'lon_2d' : lon_v,
        'lat_2d' : lat_v,
        'x'      : x,
        'y'      : y_v,
        'x_2d'   : X_v,
        'y_2d'   : Y_v,
        } 
      }

    return output
      

def find_parent_ij(outer_grid, inner_grid):
    projection = outer_grid[0]
    lat_0 = outer_grid[1]
    lon_0 = outer_grid[2]
    lat_1 = outer_grid[3]
    lat_2 = outer_grid[4]
    grid_centre_lon = outer_grid[5]
    grid_centre_lat = outer_grid[6]
    nx = outer_grid[7]
    ny = outer_grid[8]
    delta = outer_grid[9]

    from pyproj import Proj
    proj = Proj(proj=projection, lat_0=lat_0, lon_0=lon_0, lat_1=lat_1, lat_2=lat_2)
    grid_centre_x, grid_centre_y = proj(grid_centre_lon, grid_centre_lat)
    grid_x_extent = nx * delta
    grid_y_extent = ny * delta
    min_x = grid_centre_x - grid_x_extent/2.
    min_y = grid_centre_y - grid_y_extent/2.
    max_x = min_x + grid_x_extent
    max_y = min_y + grid_y_extent
    outer_x = n.arange(min_x, max_x + a_small_number, delta)
    outer_y = n.arange(min_y, max_y + a_small_number, delta)  

    projection = inner_grid[0]
    lat_0 = inner_grid[1]
    lon_0 = inner_grid[2]
    lat_1 = inner_grid[3]
    lat_2 = inner_grid[4]
    grid_centre_lon = inner_grid[5]
    grid_centre_lat = inner_grid[6]
    nx = inner_grid[7]
    ny = inner_grid[8]
    delta = inner_grid[9]
 
    grid_centre_x, grid_centre_y = proj(grid_centre_lon, grid_centre_lat)
    grid_x_extent = nx * delta
    grid_y_extent = ny * delta
    min_x = grid_centre_x - grid_x_extent/2.
    min_y = grid_centre_y - grid_y_extent/2.
    max_x = min_x + grid_x_extent
    max_y = min_y + grid_y_extent
    inner_x = n.arange(min_x, max_x + a_small_number, delta)
    inner_y = n.arange(min_y, max_y + a_small_number, delta)

    return outer_x.searchsorted(inner_x[0]), outer_y.searchsorted(inner_y[0])

def read_namelist(namelist_file):
    """read contents of namelist file and return dictionary containing all options
    
    Created 20/01/01 by Thom Chubb.

    TODO: mod_levs have a slightly different format in the namelist file, but as 
    they come last in namelist.wps I have conveniently dropped them (with a warning
    of course =) ). Whoever needs them first can come up with a fix for this.
    Untested as yet with the namelist.input file. It should work fine and may be useful 
    as a consistency check between the two files. This has been buggine me for a while.
    """

    fid=open(namelist_file)

    out_dict={}
    data = fid.readlines()
    num_lines = len(data)

    for k in range(0,num_lines):
	str = data[k].rstrip('\n').rstrip(',').split()

	if str == []:
	    pass
	elif str[0] == '':
	    pass
	elif str[0][0] == '/' or str[0][0] == '':
	    pass
	elif str[0][0] == '&':
	    # Then this line is a namelist title
	    label = str[0]

	    if label == '&mod_levs':
		print ">> WARNING: mod levels don't work yet"
		break

	    out_dict[label] ={}

	else: 
	    
	    field = str[0]
	    out_dict[label][field] = [] 

	    for k in range(2,str.__len__()):
		dat = str[k].rstrip(',')
		try:
		    dat=float(dat)
		except ValueError:
		    pass

		out_dict[label][field].append(dat) 
	    
	    # out_dict[label][field] = [] 
	    # out_dict[label][field].append(str[2:])

    return out_dict

def wrf_grid_wrapper(namelist_file='namelist.wps',nest_level=0):
    """Basic wrapper to easily visualise grids specified in namelist.wps
    
    Uses wrf.utils.read_namelist() to determine the read the appropriate variables 
    in a specified namelist file and then calls wrf.utils.wrf_grid() to define 
    the Basemap projection and show the grid over a map.

    Created 20/01/08 by Thom Chubb.

    TODO: wrf_grid() has a very basic visualisation routine and there is no 
    scope yet to view all grids simultaneously. Could use viz.utils.plot_grid() 
    with wrf.utils.wrf2latlon() as an intermediate, but the application works well
    enough for now.
    """

    # Create namelist dictionary
    nd = read_namelist(namelist_file)

    # Field editing to make python happy
    if nd['&geogrid']['map_proj'][0]=="'lambert'":
	print 'debug: modify input field lambert -> lcc' 
	nd['&geogrid']['map_proj'][0]='lcc'

    # Create wrf grid
    grd = wrf_grid(nd['&geogrid']['map_proj'][0],
			nd['&geogrid']['truelat1'][0], 
			nd['&geogrid']['truelat2'][0],
			nd['&geogrid']['stand_lon'][0],
			nd['&geogrid']['ref_lat'][0],
			nd['&geogrid']['ref_lon'][0],
			nd['&geogrid']['dx'][0],
			nd['&geogrid']['dy'][0],
			nd['&geogrid']['e_we'][nest_level],
			nd['&geogrid']['e_sn'][nest_level],
			show_mass_grid = True,
			show_stag_grids = False)


    
