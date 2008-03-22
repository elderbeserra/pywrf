"""
This function accepts the fields needed by WRF and more specifically:

based on the GFS Vtable.
Assumptions:
 - all the 3D fields need to be on isobaric levels when they are passed to the
   encoding function
 - every 3D field must be supplied with a corresponding surface field
 - the COARDS ordering of dimensions is assumed i.e. time, level, lat, lon

Usage:

"""

#############################################################################
# TODO
# - change the 'data' structure to a dictionary creating appropriate labels 
# from the dates
#############################################################################

# This script depends from a lot of external libraries:
import os, sys
import numpy as n
import pylab as p
import grib2
from string import zfill
from matplotlib.numerix import ma

# VB the following is both host and user specific hence
# VB TODO Is giving away in the source login and host names a potential
# security issue?
from socket import gethostname
hostname = gethostname()
user = os.getlogin()
if hostname == 'hn3.its.monash.edu.au':
    import PyNGL.Nio as nio
    if user == 'vbisign':
        sys.path.append('/nfs/1/home/vbisign/wrf/pywrf')   
        import viz.utils as vu
        from misc.met_utils import *
        cnvgrib_cmd = 'cnvgrib'
    elif user == 'tchubb':
        print 'Hey Thom where do you keep pywrf on this computer?'
        sys.exit()
        sys.path.append('/somewhere/pylib')
        import pywrf.viz.utils as vu
        from pywrf.misc.met_utils import *
        cnvgrib_cmd = '???'
elif hostname == 'linux450':
    # VB Sorry Thom if this is not correct ;)
    print 'Hey Thom where do you keep pywrf on this computer?'
    print 'Hey Thom how do you import pynio on this computer?'
    sys.exit()
    import PyNGL_numpy.Nio as nio
    sys.path.append('/somewhere/pylib')
    import pywrf.viz.utils as vu
    from pywrf.misc.met_utils import *
    cnvgrib_cmd = '???'
elif hostname == 'val.maths.monash.edu.au' \
    or hostname == 'valerio-bisignanesis-computer.local':
    import PyNGL_numpy.Nio as nio
    sys.path.append('/Users/val/Desktop/workspace/pywrf')
    import viz.utils as vu
    from misc.met_utils import *
    cnvgrib_cmd = '/Users/val/src/cnvgrib-1.1.3/cnvgrib'
else:
    print 'Warning: since I do not know this user/'\
      + 'hostname combination, I am not sure of ' \
      + ' where to find pywrf.viz.util, I will try\n' \
      + ' import pywrf.viz.utils as vu\n' \
      + ' from pywrf.misc.met_utils import *\n' \
      + ' import PyNGL.Nio as nio\n' \
      + " cnvgrib_cmd = '???'\n"
    import pywrf.viz.utils as vu
    from pywrf.misc.met_utils import *
    import PyNGL.Nio as nio
    cnvgrib_cmd = '???'

# VB TODO probably the data directories should be included in the above
# statement
#data_directory = '/Users/val/data/laps_nc_files/press/first_try/'
#sfc_data_directory = '/Users/val/data/laps_surface_data/'
#data_directory = '/mnt/mac/Users/val/data/laps_nc_files/press/first_try/'
#sfc_data_directory = '/mnt/mac/Users/val/data/laps_surface_data/'
#data_directory = '/media/DATA/wrf/canberra_nc_files/press/'
#sfc_data_directory = '/media/DATA/wrf/canberra_nc_files/sfc/'
#data_directory = '/media/DATA/wrf/canberra/press/'
#sfc_data_directory = '/media/DATA/wrf/canberra/sfc/'
#data_directory = '/Volumes/COMMON/wrf/canberra/press/'
#sfc_data_directory = '/Volumes/COMMON/wrf/canberra/sfc/'
#data_directory = '/Volumes/DATA/wrf/canberra/press/'
#sfc_data_directory = '/Volumes/DATA/wrf/canberra/sfc/'
#data_directory = '/Users/val/Desktop/workspace/plotting/canberra/press/'
#sfc_data_directory = '/Users/val/Desktop/workspace/plotting/canberra/sfc/'
#data_directory = \
#  '/nfs/1/home/vbisign/wrf/wps_data/pygrib2_test/canberra/press/analyses/'
#sfc_data_directory = \
#  '/nfs/1/home/vbisign/wrf/wps_data/pygrib2_test/canberra/sfc/'
data_directory = \
  '/nfs/1/home/vbisign/wrf/wps_data/run_003/press/'
sfc_data_directory = \
#sfc_data_directory = '/Users/val/Desktop/workspace/plotting/canberra/sfc/'
# It is assumed that most people will use a single nc file as the source
# of their fields. Neverthless, the field nc_file_name is available for those
# that need to use multiple data sources.

# the code is split in two parts... the first part gathers all the
# information that is specific to your data source... stuff like the location
# of the netcdf files, the name of the variables of interest in your files,
# and any special treatment of the fields like changing units, scaling the
# data or anything that is needed to get the data in the form required by the
# encoding function.
# This translates to simply having to generate a fairly pedantic dictionary
# that spells out what needs to be done with the fields during the grib
# encoding
# the second part is the encoding function itself in which it is assumed the
# data is passed in a dictionary containing the data as numpy arrays in the 
# right units

# Unfortunately it proves hard to both generalize the data structure and
# maintain flexibility... hence at the price of making the code longer I am
# going to process one field at the time so that it is clear to anyone what
# the steps that are (may be) involved are...
# I still believe it is a good idea to collect all the information needed in
# one big dictionary for ease of retrieval


#############################################################################
# Preliminaries
# This is the data structure that I will pass to the encoder
data = []

# to generate a grib2 file using Jeffrey Whitaker's grib2 class, it is necessary
# to provide the constructor with two pieces of information:
# - the discipline of the data to be encoded
# - the identification section
# this means the class will need to be reconstructed every time the discipline
# changes


def missing(octects):
    # To say that something is not available in a grib field they use a 
    # sequence of 1s as long as the space allocated (number of octects).
    # When this sequence is read as an integer, it corresponds to the maximum
    # representable number given the octects*8 bits of memory, hence:
    # PS: refer to the grib documentation to know how many octects are
    # allocated to any particular field
    return 2**(8*octects) - 1

def prepare_sfc_data(f, fv):
    # work out which file to pick based on the base date of the first file
    # this is based on the (fragile) assumption that nobody changed the file
    # names from the standard laps...
    # this was needed for the first batch of data...
    valid_date = fv['valid_date'].get_value()[0]
    valid_time = fv['valid_time'].get_value()[0]
    # the following is old stuff
    if 0:
        # the data in each file 'starts' at 12:00pm (UTC)
        if valid_time < 1200:
            valid_date -= 1
        #sfc_file = 'laps' + str(valid_date)[4:] + '_surface.nc'
        print sfc_data_directory + sfc_file
    # this is new and (hopefully) better
    if 1:
        # there is one sfc data file for each of the analysis containing the
        # full 73 hours of diagnostic fields
        # in the simplified case that we do not span analysis -> we start again
        # with a different file at every analysis:
        base_date = fv['base_date'].get_value()
        base_time = fv['base_time'].get_value()
        sfc_file = 'regprog.laps_pt375.' \
          + str(base_date)[4:] + '.' + str(base_time).zfill(4)[:2] \
          + '.slvfld.n3.nc'
    # if we use data from a single forecast we just name the appropriate file
    # sfc_file = 'laps0117_surface.nc'
    sfc_f = nio.open_file(sfc_data_directory + sfc_file)
    sfc_fv = sfc_f.variables
    # I need this to work out which time to pick in the file...
    sfc_valid_time = sfc_fv['valid_time'].get_value()
    time_idx = sfc_valid_time.tolist().index(valid_time)
    return sfc_file, sfc_f, sfc_fv, time_idx

#############################################################################

#############################################################################
# this variables remain unchanged across the different times
# BTW these (with the exception of the significance_ref_time) should not
# affect WRF

# Id of originating centre (Common Code Table C-1) -> 1 for Melbourne
# Id of originating centre (Common Code Table C-1) -> 7 for NCEP
centre = 7 
# Id of orginating sub-centre (local table) -> 0 for no subcentre
# Id of orginating sub-centre (local table) -> 3 for NCEP Central Operations
subcentre = 0
# GRIB Master Tables Version Number (Code Table 1.0)
master_table = 1
# GRIB Local Tables Version Number (Code Table 1.1) -> 0 for not used
# GRIB Local Tables Version Number (Code Table 1.1) -> 1 for let's see what happens
local_table = 1
# Significance of Reference Time (Code Table 1.2) -> 1 for start of forecast
significance_ref_time = 1
# Production status of data (Code Table 1.3) -> 2 for research product
production_status = 2
# idsect[12]=Type of processed data (Code Table 1.4) -> 1 -> forecast product
type_processed_data = 1

#############################################################################

# the first dimension of data will be the times to be processed
# you will have to modify this to suit your data
# this could mean just having a time index to get different times from the
# same file or just pick different files for different times.
files = os.listdir(data_directory)
#for file in files[:2]:
for file in files:
    f = nio.open_file(data_directory + file)
    fv = f.variables
    base_date = str(fv['base_date'].get_value())
    year = int(base_date[0:4])
    month = int(base_date[4:6])
    day = int(base_date[6:8])
    # padding needed to use a string method on an integer
    base_time = zfill(str(fv['base_time'].get_value()),4)
    hour = int(base_time[0:2])
    minute =  int(base_time[2:])
    second = 0
    # I did not generate this in the preamble above because I needed to wait
    # for the base date and time info from the netcdf files...
    # NB the order of the variables in the list MATTERS!
    idsect = [
      centre,
      subcentre,
      master_table,
      local_table,
      significance_ref_time,
      year,
      month,
      day,
      hour,
      minute,
      second,
      production_status,
      type_processed_data
    ]

    print file
    # gdsinfo - Sequence containing information needed for the grid definition 
    #     section.
    gdsinfo = n.zeros((5,), dtype=int)
    # gdsinfo[0] = Source of grid definition (see Code Table 3.0) -> 0 -> will be
    #     specified in octects 13-14
    gdsinfo[0] = 0
    # gdsinfo[1] = Number of grid points in the defined grid.
    #     TODO make it more general
    #lat = fv['lat'].get_value()
    # the idx was chosen to stop the data at 55S and avoid the rubbish far south
    lat_cheat_idx = 30
    # 1 for the full range
    lon_cheat_idx = 1
    lat = fv['lat'].get_value()[lat_cheat_idx:]
    lon = fv['lon'].get_value()[:-lon_cheat_idx]
    gdsinfo[1] = len(lat)*len(lon)
    # gdsinfo[2] = Number of octets needed for each additional grid points defn. 
    #     Used to define number of points in each row for non-reg grids 
    #     (=0 for regular grid).
    gdsinfo[2] = 0
    # gdsinfo[3] = Interp. of list for optional points defn (Code Table 3.11)
    #     0 -> There is no appended list
    gdsinfo[3] = 0
    # gdsinfo[4] = Grid Definition Template Number (Code Table 3.1)
    #     0 -> Latitude/Longitude (or equidistant cylindrical, or Plate Carree)
    gdsinfo[4] = 0
    
    # gdtmpl - Contains the data values for the specified Grid Definition 
    #     Template ( NN=gds[4] ). Each element of this integer array contains 
    #     an entry (in the order specified) of Grid Definition Template 3.NN
    gdtmpl = n.zeros((19,), dtype=n.int64)
    # In the following we refer to the lat/lon template (grid template 3.0)
    # Oct: 15      Shape of the Earth (See Table 3.2)
    #     0 -> spherical with radius = 6,367,470.0 m
    #     TODO double check the radius they used in Tan's code
    gdtmpl[0] = 0
    # Oct: 16      Scale Factor of radius of spherical Earth
    #     TODO hopefully (given the last option) the next 5 are ignored...
    #     if 1 octet is the basic unit then a value of 255 should be equivalent
    #     to all bits set to 1 -> this could be broken if more than one octect
    #     is assigned -> maybe I should use the maximum represantable number
    #     tailored to the number of octects?
    #     are there any rules that make the encoder ignore automatically certain
    #     options e.g. the oblate spheroid one when a spherical shape is assumed
    gdtmpl[1] = missing(1)
    # Oct: 17-20   Scale value of radius of spherical Earth
    # NB they (grib2 std) interprets all bits set to 1 as a missing value
    #     hence the value is base*(word_width*words_assigned) - 1
    #     base is 2 since we are working in binary
    #     word_width is 8 since they use octects
    #     words_assigned is the number of words reserved to represent the value
    #     -1 gives us the greatest representable number since we start from 0
    gdtmpl[2] = missing(1)
    # Oct: 21      Scale factor of major axis of oblate spheroid Earth
    gdtmpl[3] = missing(1)
    # Oct: 22-25   Scaled value of major axis of oblate spheroid Earth
    gdtmpl[4] = missing(4)
    # Oct: 26      Scale factor of minor axis of oblate spheroid Earth
    gdtmpl[5] = missing(1)
    # Oct: 27-30   Scaled value of minor axis of oblate spheroid Earth
    gdtmpl[6] = missing(4)
    # Oct: 31-34   Number of points along a parallel     
    gdtmpl[7] = len(lon)
    # Oct: 35-38   Number of points along a meridian
    gdtmpl[8] = len(lat)
    # Oct: 39-42   Basic angle of the initial production domain (see Note 1)
    gdtmpl[9] = missing(4)
    # Oct: 43-46   Subdivisions of basic angle used to define extreme longitudes 
    #             and latitudes, and direction increments (see Note 1)
    gdtmpl[10] = missing(4)
    # Oct: 47-50   La1 - Latitude of first grid point (see Note 1)        
    #gdtmpl.append(in_vars['lat'].get_value()[-1]*1e6)
    gdtmpl[11] = lat[-1]*1e6
    #gdtmpl[11] = lat[0]*1e6
    # Oct: 51-54   Lo1 - Longitude of first grid point (see Note 1)
    #gdtmpl.append(in_vars['lon'].get_value()[0]*1e6)
    gdtmpl[12] = lon[0]*1e6
    # Oct: 55      Resolution and component flags (see Table 3.3)
    #    This should set all the bits to 0 meaning 
    #    bit        value   meaning
    #    1-2        ?       Reserved
    #    3          0       i direction increments not given
    #    4	        0       j direction increments not given
    #    5          0       Resolved u and v components of vector quantities 
    #                       relative to easterly and northerly directions
    #    6-8        0       Reserved - set to zero
    #gdtmpl[13] = 0 
    #gdtmpl[13] = 12 
    gdtmpl[13] = 48 
    # Oct: 56-59   La2 - Latitude of last grid point (see Note1)     
    #gdtmpl.append(in_vars['lat'].get_value()[0]*1e6)
    gdtmpl[14] = lat[0]*1e6
    #gdtmpl[14] = lat[-1]*1e6
    # Oct: 60-63   Lo2 - Longitude of last grid point (see Note 1)
    #gdtmpl.append(in_vars['lon'].get_value()[-1]*1e6)
    gdtmpl[15] = lon[-1]*1e6
    # Oct: 64-67   i (west-east) Direction increment (see Note1)     
    #gdtmpl[16] = missing(4)
    gdtmpl[16] = 0.375*1e6
    # Oct: 68-71   j (south-north) Direction increment (see Note1)
    #gdtmpl[17] = missing(4)
    gdtmpl[17] = 0.375*1e6
    #gdtmpl[17] = -0.375*1e6
    #gdtmpl[17] = 33143
    # Oct: 72      Scanning mode (flags see Table 3.4)
    #    TODO double check if this is legit considering the flags 
    #    This should set all the bits to 0 meaning 
    #    bit        value   meaning
    #    1          0       Points in the first row or column scan in +i (+x) 
    #                       direction
    #    2          0       Points in the first row or column scan in -i (-x) 
    #                       direction
    #    3          0       Adjacent pints in i(x) direction are consecutive
    #    4	        0       All rows scan in the same direction
    #    5-8        0       Reserved (set to zero)
    #    NB I suspect the 3 bit might be important for the 'majoring' of the array?
    #           	# 2^7	2^6	2^5	2^4	2^3	2^2	2^1	2^0
    #           	# 1	2	3	4	5	6	7	8
    gdtmpl[18] = 0 	# 0	0	0	0	0	0	0	0
    #gdtmpl[18] = 16 	# 0	0	0	1	0	0	0	0
    #gdtmpl[18] = 32 	# 0	0	1	0	0	0	0	0
    #gdtmpl[18] = 64 	# 0	1	0	0	0	0	0	0
    #gdtmpl[18] = 128 	# 1	0	0	0	0	0	0	0

    #gdtmpl[18] = 192 	# 1	1	0	0	0	0	0	0
    #gdtmpl[18] = 160 	# 1	0	1	0	0	0	0	0
    #gdtmpl[18] = 144 	# 1	0	0	1	0	0	0	0
    #gdtmpl[18] = 96 	# 0	1	1	0	0	0	0	0
    #gdtmpl[18] = 80 	# 0	1	0	1	0	0	0	0
    #gdtmpl[18] = 48  	# 0	0	1	1	0	0	0	0

    #gdtmpl[18] = 224  	# 1	1	1	0	0	0	0	0
    #gdtmpl[18] = 208  	# 1	1	0	1	0	0	0	0
    #gdtmpl[18] = 176  	# 1	0	1	1	0	0	0	0
    #gdtmpl[18] = 112  	# 0	1	1	1	0	0	0	0

    #gdtmpl[18] = 240  	# 1	1	1	1	0	0	0	0

    data.append(
      {
        'idsect' : idsect,
        'gdsinfo' : gdsinfo,
        'gdtmpl' : gdtmpl,  
        'discipline':
          {
            # the first number in the list is the code for the discipline
            'atmos': 
               {
                 'code':0,
                 'fields':{'2d':{}, '3d':{}}
               },
            'ocean': 
               {
                 'code':10,
                 'fields':{'2d':{}, '3d':{}}
               },
            'land': 
               {
                 'code':2,
                 'fields':{'2d':{}, '3d':{}}
               }
          }
      }
    )
    # this is the part where we extract the data from the nc files and put
    # into the structure we will pass later to the encoder
    # let's start with the atmospheric products specifically the surface
    # fields
    # each variables comes with all the metadata necessary to do its grib
    # encoding using the addfield method of the encoder class i.e. pdtnum,
    # pdtmpl, drtnum, drtmpl, field

    # -1 stands for the last one (time) added
    dummy = data[-1]['discipline']['atmos']['fields']['2d']
    
    # Let's prepare the surface pressure data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 3 -> mass
    pdtmpl[0] = 3
    # Oct 11: Parameter number (see Code table 4.2). 0 -> pressure
    pdtmpl[1] = 0
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 1 -> ground or water surface
    pdtmpl[9] = 1
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    pdtmpl[11] = 0
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)

    # drtnum - Data Representation Template Number (see Code Table 5.0).
    #          0 -> grid point data - simple packing
    drtnum = 0
    # drtmpl - Sequence with the data values for the specified Data Representation 
    #     Template (N=drtnum). Each element of this integer array contains an 
    #     entry (in the order specified) of Data Representation Template 5.N Note 
    #     that some values in this template (eg. reference values, number of 
    #     bits, etc...) may be changed by the data packing algorithms. Use this to 
    #     specify scaling factors and order of spatial differencing, if desired.
    drtmpl = n.zeros((5,),dtype=int)
    # NB to make sense of the following remember (from the wmo std):
    # The original data value Y (in the units of code table 4.2) can be recovered 
    # with the formula:
    # Y * 10^D= R + (X1+X2) * 2^E
    # For simple packing and all spectral data
    #     E = Binary scale factor,
    #     D = Decimal scale factor
    #     R = Reference value of the whole field,
    #     X1 = 0,
    #     X2 = Scaled (encoded) value.
    
    # Oct: 12-15. Reference value (R) (IEEE 32-bit floating-point value)
    drtmpl[0] = 0
    # Oct: 16-17. Binary scale factor (E)
    drtmpl[1] = 0
    # Oct: 18-19. Decimal scale factor (D)
    drtmpl[2] = 0
    # Oct: 20. Number of bits used for each packed value for simple packing, or 
    drtmpl[3] = 24
    # for each group reference value for complex packing or spatial differencing
    # Oct: 21. Type of original field values (see Code Table 5.1)
    # 0 -> floating point
    drtmpl[4] = 0

    cheat = fv['sfc_pres'].get_value()[0,lat_cheat_idx:,:-lon_cheat_idx]
    field = n.zeros(cheat.shape)
    for lat_idx in range(cheat.shape[0]):
        field[-(lat_idx+1)] = cheat[lat_idx]

    if 0:
        dummy['sfc_pres'] = {
            'long_name' : 'surface pressure',
            'pdtnum' : pdtnum, 
            'pdtmpl' : pdtmpl, 
            'drtnum' : drtnum,
            'drtmpl' : drtmpl,
            #'field' : fv['sfc_pres'].get_value()[0]
            #'field' : fv['sfc_pres'].get_value()[0,lat_cheat_idx:,:-lon_cheat_idx]
            'field' : field
          }

    # Let's prepare the mean sea level pressure data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 3 -> mass
    pdtmpl[0] = 3
    # Oct 11: Parameter number (see Code table 4.2). 1 -> MSLP
    pdtmpl[1] = 1
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 101 -> sea level
    pdtmpl[9] = 101
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    pdtmpl[11] = 0
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)
 
    # NB I am using the same drtnum and drtmpl for all data so I want
    # replicate it in the following fields if you need to change this on a
    # per-variable basis this is the place to redefine it. If you choose to do
    # so then make sure you define one for each variable that follows or be
    # aware that the last one defined will apply to all that follows.

    cheat = fv['mslp'].get_value()[0,lat_cheat_idx:,:-lon_cheat_idx]
    # to go from mb to Pa
    #cheat *= 100.
    cheat = cheat * 100.
    field = n.zeros(cheat.shape)
    for lat_idx in range(cheat.shape[0]):
        field[-(lat_idx+1)] = cheat[lat_idx]

    dummy['mslp'] = {
        'long_name' : 'mean sea level pressure',
        # in my data it comes in mb and I want in  Pascals
        'pdtnum' : pdtnum, 
        'pdtmpl' : pdtmpl, 
        'drtnum' : drtnum,
        'drtmpl' : drtmpl,
        #'field' : fv['mslp'].get_value()[0]*100.
        #'field' : fv['mslp'].get_value()[0,lat_cheat_idx:,:-lon_cheat_idx]*100.
        'field' : field
      }

    # Let's prepare the 2m temperature data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 0 -> temperature
    pdtmpl[0] = 0
    # Oct 11: Parameter number (see Code table 4.2). 0 -> temperature
    pdtmpl[1] = 0
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 1 -> ground or water surface
    pdtmpl[9] = 103
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    # -> 2m above ground
    pdtmpl[11] = 2
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)
 
    # NB I am using the same drtnum and drtmpl for all data so I want
    # replicate it in the following fields if you need to change this on a
    # per-variable basis this is the place to redefine it. If you choose to do
    # so then make sure you define one for each variable that follows or be
    # aware that the last one defined will apply to all that follows.

    cheat = fv['sfc_temp'].get_value()[0,lat_cheat_idx:,:-lon_cheat_idx]
    field = n.zeros(cheat.shape)
    for lat_idx in range(cheat.shape[0]):
        field[-(lat_idx+1)] = cheat[lat_idx]

    dummy['sfc_temp'] = {
        'long_name' : 'surface temperature',
        'pdtnum' : pdtnum, 
        'pdtmpl' : pdtmpl, 
        'drtnum' : drtnum,
        'drtmpl' : drtmpl,
        #'field' : fv['sfc_temp'].get_value()[0]
        #'field' : fv['sfc_temp'].get_value()[0,lat_cheat_idx:,:-lon_cheat_idx]
        'field' : field
      }

    # Let's prepare the 2m temperature data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 0 -> temperature
    pdtmpl[0] = 0
    # Oct 11: Parameter number (see Code table 4.2). 0 -> temperature
    pdtmpl[1] = 0
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 1 -> surface
    pdtmpl[9] = 1
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    # -> 2m above ground
    pdtmpl[11] = 0
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)
 
    # NB I am using the same drtnum and drtmpl for all data so I want
    # replicate it in the following fields if you need to change this on a
    # per-variable basis this is the place to redefine it. If you choose to do
    # so then make sure you define one for each variable that follows or be
    # aware that the last one defined will apply to all that follows.

    cheat = fv['skin_temp'].get_value()[0,lat_cheat_idx:,:-lon_cheat_idx]
    field = n.zeros(cheat.shape)
    for lat_idx in range(cheat.shape[0]):
        field[-(lat_idx+1)] = cheat[lat_idx]

    dummy['skin_temp'] = {
        'long_name' : 'skin temperature',
        'pdtnum' : pdtnum, 
        'pdtmpl' : pdtmpl, 
        'drtnum' : drtnum,
        'drtmpl' : drtmpl,
        #'field' : fv['skin_temp'].get_value()[0]
        #'field' : fv['skin_temp'].get_value()[0,lat_cheat_idx:,:-lon_cheat_idx]
        'field' : field
      }

    # the next fields will come from a different file...
    sfc_file, sfc_f, sfc_fv, time_idx = prepare_sfc_data(f,fv)
    # DEBUG
    #p.figure()
    #p.contour(dummy['skin_temp']['field'])
    #p.figure()
    #p.contour(sfc_fv['skin_temp'].get_value()[time_idx])

    # Let's prepare the 2m RH data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 1 -> moisture
    pdtmpl[0] = 1
    # Oct 11: Parameter number (see Code table 4.2). 1 -> RH
    pdtmpl[1] = 1
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 103 -> fixed height above ground
    pdtmpl[9] = 103
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    # -> 2m above ground
    pdtmpl[11] = 2
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)
 
    # NB I am using the same drtnum and drtmpl for all data so I want
    # replicate it in the following fields if you need to change this on a
    # per-variable basis this is the place to redefine it. If you choose to do
    # so then make sure you define one for each variable that follows or be
    # aware that the last one defined will apply to all that follows.

    # The surface files do not come with mixing ratio, but we can calculate
    # that from the dew point temperature...
    #dew_point = sfc_fv['scrn_dewpt'].get_value()[time_idx]
    dew_point = sfc_fv['scrn_dewpt'].get_value()[time_idx,lat_cheat_idx:,:-lon_cheat_idx]
    vapour_pressure = goff_gratch(dew_point)
    #sfc_temp = sfc_fv['scrn_temp'].get_value()[time_idx]
    sfc_temp = sfc_fv['scrn_temp'].get_value()[time_idx,lat_cheat_idx:,:-lon_cheat_idx]
    sat_vapour_pressure = goff_gratch(sfc_temp)
    #sat_vapour_pressure = ma.masked_where(sat_vapour_pressure == 0., 
    #  sat_vapour_pressure)
    sfc_rh = vapour_pressure / sat_vapour_pressure * 100
    del sfc_temp, vapour_pressure, dew_point, sat_vapour_pressure

    cheat = sfc_rh
    field = n.zeros(cheat.shape)
    for lat_idx in range(cheat.shape[0]):
        field[-(lat_idx+1)] = cheat[lat_idx]

    dummy['sfc_rh'] = {
        'long_name' : '2m relative humidity',
        'pdtnum' : pdtnum, 
        'pdtmpl' : pdtmpl, 
        'drtnum' : drtnum,
        'drtmpl' : drtmpl,
        #'field' : sfc_rh.filled(fill_value=0.)
        #'field' : sfc_rh
        'field' : field
      }

    # Let's prepare the zonal 10m wind data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 2 -> momentum
    pdtmpl[0] = 2
    # Oct 11: Parameter number (see Code table 4.2). 2 -> u-component of wind
    pdtmpl[1] = 2
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 103 -> fixed height above ground
    pdtmpl[9] = 103
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    # -> 10m above ground
    pdtmpl[11] = 10
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)
 
    # NB I am using the same drtnum and drtmpl for all data so I want
    # replicate it in the following fields if you need to change this on a
    # per-variable basis this is the place to redefine it. If you choose to do
    # so then make sure you define one for each variable that follows or be
    # aware that the last one defined will apply to all that follows.

    cheat = sfc_fv['q10m_zonal_wnd'].get_value()[time_idx,lat_cheat_idx:,:-lon_cheat_idx]
    field = n.zeros(cheat.shape)
    for lat_idx in range(cheat.shape[0]):
        field[-(lat_idx+1)] = cheat[lat_idx]

    dummy['sfc_zonal_wnd'] = {
        'long_name' : '10m zonal wind',
        'pdtnum' : pdtnum, 
        'pdtmpl' : pdtmpl, 
        'drtnum' : drtnum,
        'drtmpl' : drtmpl,
        #'field' : sfc_fv['q10m_zonal_wnd'].get_value()[time_idx]
        #'field' : sfc_fv['q10m_zonal_wnd'].get_value()[time_idx,lat_cheat_idx:,:-lon_cheat_idx]
        'field' : field
      }

    # Let's prepare the meridional 10m wind data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 2 -> momentum
    pdtmpl[0] = 2
    # Oct 11: Parameter number (see Code table 4.2). 3 -> v-component of wind
    pdtmpl[1] = 3
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 103 -> fixed height above ground
    pdtmpl[9] = 103
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    # -> 10m above ground
    pdtmpl[11] = 10
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)
 
    # NB I am using the same drtnum and drtmpl for all data so I want
    # replicate it in the following fields if you need to change this on a
    # per-variable basis this is the place to redefine it. If you choose to do
    # so then make sure you define one for each variable that follows or be
    # aware that the last one defined will apply to all that follows.

    cheat = sfc_fv['q10m_merid_wnd'].get_value()[time_idx,lat_cheat_idx:,:-lon_cheat_idx]
    field = n.zeros(cheat.shape)
    for lat_idx in range(cheat.shape[0]):
        field[-(lat_idx+1)] = cheat[lat_idx]

    dummy['sfc_merid_wnd'] = {
        'long_name' : '10m meridional wind',
        'pdtnum' : pdtnum, 
        'pdtmpl' : pdtmpl, 
        'drtnum' : drtnum,
        'drtmpl' : drtmpl,
        #'field' : sfc_fv['q10m_merid_wnd'].get_value()[time_idx]
        #'field' : sfc_fv['q10m_merid_wnd'].get_value()[time_idx,lat_cheat_idx:,:-lon_cheat_idx]
        'field' : field
      }
  
    if 1:
        # Let's prepare the topography data and metadata
        # Product Definition Template Number (see Code Table 4.0)
        #     0 -> Analysis or forecast at a horizontal level or in a 
        #          horizontal layer at a point in time.  (see Template 4.0)
        pdtnum = 0
        # array of zeros (int) of size 15
        pdtmpl = n.zeros((15,),dtype=int)
        # Follows the definition of template 4.0 
        # Oct 10: Parameter category (see Code table 4.1). 3 -> mass
        pdtmpl[0] = 3
        # Oct 11: Parameter number (see Code table 4.2). 5 -> geopotential height
        pdtmpl[1] = 5
        # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
        pdtmpl[2] = 2
        # Oct 13: Background generating process identifier 
        # (defined by originating centre)
        pdtmpl[3] = missing(1)
        # Oct 14: Analysis or forecast generating process identified 
        # (defined by originating centre)
        pdtmpl[4] = missing(1)
        # Oct 15-16: Hours of observational data cutoff after reference time 
        pdtmpl[5] = missing(2)
        # Oct 17: Minutes of observational data cutoff after reference time 
        pdtmpl[6] = missing(1)
        # Oct 18: Indicator of unit of time range (see Code table 4.4)
        # 1 -> hours
        # NB WPS expects this to be hours!
        pdtmpl[7] = 1
        # Oct 19-22: Forecast time in units defined by octet 18
        pdtmpl[8] = fv['forc_hrs'].get_value()[0]
        # Oct 23: Type of first fixed surface (see Code table 4.5)
        # 1 -> surface
        pdtmpl[9] = 1
        # Oct 24: Scale factor of first fixed surface
        pdtmpl[10] = 0
        # Oct 25-28: Scaled value of first fixed surface
        pdtmpl[11] = 0
        # Oct 29: Type of second fixed surface (see Code table 4.5)
        pdtmpl[12] = missing(1)
        # Oct 30: Scale factor of second fixed surface
        pdtmpl[13] = missing(1)
        # Oct 31-34: Scaled value of second fixed surfaces
        pdtmpl[14] = missing(3)
     
        # NB I am using the same drtnum and drtmpl for all data so I want
        # replicate it in the following fields if you need to change this on a
        # per-variable basis this is the place to redefine it. If you choose to do
        # so then make sure you define one for each variable that follows or be
        # aware that the last one defined will apply to all that follows.
    
        cheat = fv['topog'].get_value()[lat_cheat_idx:,:-lon_cheat_idx]
        field = n.zeros(cheat.shape)
        for lat_idx in range(cheat.shape[0]):
            field[-(lat_idx+1)] = cheat[lat_idx]

        dummy['topog'] = {
            'long_name' : 'topography',
            'pdtnum' : pdtnum, 
            'pdtmpl' : pdtmpl, 
            'drtnum' : drtnum,
            'drtmpl' : drtmpl,
            #'field' : fv['topog'].get_value()
            #'field' : fv['topog'].get_value()[lat_cheat_idx:,:-lon_cheat_idx]
            'field' : field
          }

    if 1:
        # Let's prepare the water equivalent accumulated snow depth
        # data and metadata
        # Product Definition Template Number (see Code Table 4.0)
        #     0 -> Analysis or forecast at a horizontal level or in a 
        #          horizontal layer at a point in time.  (see Template 4.0)
        pdtnum = 0
        # array of zeros (int) of size 15
        pdtmpl = n.zeros((15,),dtype=int)
        # Follows the definition of template 4.0 
        # Oct 10: Parameter category (see Code table 4.1). 1 -> moisture
        pdtmpl[0] = 1
        # Oct 11: Parameter number (see Code table 4.2). 13 -> 
        #   water equivalent of accumulated snow depth
        pdtmpl[1] = 13
        # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
        pdtmpl[2] = 2
        # Oct 13: Background generating process identifier 
        # (defined by originating centre)
        pdtmpl[3] = missing(1)
        # Oct 14: Analysis or forecast generating process identified 
        # (defined by originating centre)
        pdtmpl[4] = missing(1)
        # Oct 15-16: Hours of observational data cutoff after reference time 
        pdtmpl[5] = missing(2)
        # Oct 17: Minutes of observational data cutoff after reference time 
        pdtmpl[6] = missing(1)
        # Oct 18: Indicator of unit of time range (see Code table 4.4)
        # 1 -> hours
        # NB WPS expects this to be hours!
        pdtmpl[7] = 1
        # Oct 19-22: Forecast time in units defined by octet 18
        pdtmpl[8] = fv['forc_hrs'].get_value()[0]
        # Oct 23: Type of first fixed surface (see Code table 4.5)
        # 1 -> surface
        pdtmpl[9] = 1
        # Oct 24: Scale factor of first fixed surface
        pdtmpl[10] = 0
        # Oct 25-28: Scaled value of first fixed surface
        pdtmpl[11] = 0
        # Oct 29: Type of second fixed surface (see Code table 4.5)
        pdtmpl[12] = missing(1)
        # Oct 30: Scale factor of second fixed surface
        pdtmpl[13] = missing(1)
        # Oct 31-34: Scaled value of second fixed surfaces
        pdtmpl[14] = missing(3)
     
        # NB I am using the same drtnum and drtmpl for all data so I want
        # replicate it in the following fields if you need to change this on a
        # per-variable basis this is the place to redefine it. If you choose to do
        # so then make sure you define one for each variable that follows or be
        # aware that the last one defined will apply to all that follows.
    
        #dummy_field = n.zeros((len(lon),len(lat)),dtype=float)
        dummy_field = n.zeros((len(lat),len(lon)),dtype=float)
        dummy['weasd'] = {
            'long_name' : 'water equivalent of accumulated snow depth',
            'pdtnum' : pdtnum, 
            'pdtmpl' : pdtmpl, 
            'drtnum' : drtnum,
            'drtmpl' : drtmpl,
            'field' : dummy_field
          }
       
    #######################################################################
    # Finished with the 2d fields
    # Let's move on to the 3d ones
    #######################################################################
    
    dummy = data[-1]['discipline']['atmos']['fields']['3d']
    
    # Let's prepare the air temperature data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 0 -> temperature
    pdtmpl[0] = 0
    # Oct 11: Parameter number (see Code table 4.2). 0 -> temperature
    pdtmpl[1] = 0
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 1 -> isobaric level
    pdtmpl[9] = 100
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    # -> missing as there is will be defined later when the data is sliced in
    # horizontal slabs for its encoding
    pdtmpl[11] = missing(3)
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)
 
    # NB I am using the same drtnum and drtmpl for all data so I want
    # replicate it in the following fields if you need to change this on a
    # per-variable basis this is the place to redefine it. If you choose to do
    # so then make sure you define one for each variable that follows or be
    # aware that the last one defined will apply to all that follows.

    cheat = fv['air_temp'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx]
    field = n.zeros(cheat.shape)
    for lvl_idx in range(cheat.shape[0]):
        for lat_idx in range(cheat.shape[1]):
            field[lvl_idx,-(lat_idx+1)] = cheat[lvl_idx,lat_idx]

    dummy['air_temp'] = {
        'long_name' : '3D air temperature',
        'pdtnum' : pdtnum, 
        'pdtmpl' : pdtmpl, 
        'drtnum' : drtnum,
        'drtmpl' : drtmpl,
        #'field' : fv['air_temp'].get_value()[0],
        #'field' : fv['air_temp'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx],
        'field' : field,
        #'field' : '',
        'lvl' : fv['lvl'].get_value()
      }

    # Let's prepare the 3d relative humidity data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 1 -> moisture
    pdtmpl[0] = 1
    # Oct 11: Parameter number (see Code table 4.2). 1 -> RH
    pdtmpl[1] = 1
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 1 -> isobaric level
    pdtmpl[9] = 100
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    # -> missing as there is will be defined later when the data is sliced in
    # horizontal slabs for its encoding
    pdtmpl[11] = missing(3)
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)
 
    # NB I am using the same drtnum and drtmpl for all data so I want
    # replicate it in the following fields if you need to change this on a
    # per-variable basis this is the place to redefine it. If you choose to do
    # so then make sure you define one for each variable that follows or be
    # aware that the last one defined will apply to all that follows.

    #mix_rto = fv['mix_rto'].get_value()[0]
    mix_rto = fv['mix_rto'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx]
    spec_hum = n.zeros(mix_rto.shape,dtype=float)
    spec_hum = mix_rto / (mix_rto + 1)
    #temp = fv['air_temp'].get_value()[0]
    temp = fv['air_temp'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx]
    #print temp.min(), temp.max()
    #print temp.min()-273.16, temp.max()-273.16
    sat_vapour_pres = goff_gratch(temp)
    #print sat_vapour_pres.min(), sat_vapour_pres.max()
    pres = fv['lvl'].get_value()
    sat_spec_hum = n.zeros(mix_rto.shape,dtype=float)
    for lvl_idx in range(len(pres)):
        sat_spec_hum[lvl_idx] = 0.622 * sat_vapour_pres[lvl_idx] \
          / pres[lvl_idx]
    #sat_spec_hum = ma.masked_where(sat_spec_hum <= 1.e-9,
    #sat_spec_hum = ma.masked_where(sat_spec_hum <= 1.e-6,
      #sat_spec_hum)
    rh = n.zeros(mix_rto.shape,dtype=float)
    rh = spec_hum / sat_spec_hum * 100
    #sys.exit()
    del mix_rto, spec_hum, temp, sat_vapour_pres, pres, sat_spec_hum
    
    #cheat = rh
    #field = n.zeros(cheat.shape)
    #for lat_idx in range(cheat.shape[0]):
    #    field[-(lat_idx+1)] = cheat[lat_idx]

    cheat = rh
    field = n.zeros(cheat.shape)
    for lvl_idx in range(cheat.shape[0]):
        for lat_idx in range(cheat.shape[1]):
            field[lvl_idx,-(lat_idx+1)] = cheat[lvl_idx,lat_idx]

    dummy['rh'] = {
        'long_name' : '3D relative humidity',
        'pdtnum' : pdtnum, 
        'pdtmpl' : pdtmpl, 
        'drtnum' : drtnum,
        'drtmpl' : drtmpl,
        #'field' : rh.filled(fill_value=0.),
        #'field' : rh,
        'field' : field,
        'lvl' : fv['lvl'].get_value()
      }
    #sys.exit()

    # Let's prepare the 3d zonal wind data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 2 -> momentum
    pdtmpl[0] = 2
    # Oct 11: Parameter number (see Code table 4.2). 2 -> u-component of the
    # wind
    pdtmpl[1] = 2
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 1 -> isobaric level
    pdtmpl[9] = 100
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    # -> missing as there is will be defined later when the data is sliced in
    # horizontal slabs for its encoding
    pdtmpl[11] = missing(3)
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)
 
    # NB I am using the same drtnum and drtmpl for all data so I want
    # replicate it in the following fields if you need to change this on a
    # per-variable basis this is the place to redefine it. If you choose to do
    # so then make sure you define one for each variable that follows or be
    # aware that the last one defined will apply to all that follows.

    cheat = fv['zonal_wnd'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx]
    field = n.zeros(cheat.shape)
    for lvl_idx in range(cheat.shape[0]):
        for lat_idx in range(cheat.shape[1]):
            field[lvl_idx,-(lat_idx+1)] = cheat[lvl_idx,lat_idx]

    dummy['zonal_wnd'] = {
        'long_name' : '3D zonal wind',
        'pdtnum' : pdtnum, 
        'pdtmpl' : pdtmpl, 
        'drtnum' : drtnum,
        'drtmpl' : drtmpl,
        #'field' : fv['zonal_wnd'].get_value()[0],
        #'field' : fv['zonal_wnd'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx],
        'field' : field,
        #'field' : '',
        'lvl' : fv['lvl'].get_value()
      }

    # Let's prepare the 3d zonal wind data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 2 -> momentum
    pdtmpl[0] = 2
    # Oct 11: Parameter number (see Code table 4.2). 2 -> v-component of the
    # wind
    pdtmpl[1] = 3
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 1 -> isobaric level
    pdtmpl[9] = 100
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    # -> missing as there is will be defined later when the data is sliced in
    # horizontal slabs for its encoding
    pdtmpl[11] = missing(3)
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)
 
    # NB I am using the same drtnum and drtmpl for all data so I want
    # replicate it in the following fields if you need to change this on a
    # per-variable basis this is the place to redefine it. If you choose to do
    # so then make sure you define one for each variable that follows or be
    # aware that the last one defined will apply to all that follows.

    cheat =fv['merid_wnd'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx] 
    field = n.zeros(cheat.shape)
    for lvl_idx in range(cheat.shape[0]):
        for lat_idx in range(cheat.shape[1]):
            field[lvl_idx,-(lat_idx+1)] = cheat[lvl_idx,lat_idx]

    dummy['merid_wnd'] = {
        'long_name' : '3D meridional wind',
        'pdtnum' : pdtnum, 
        'pdtmpl' : pdtmpl, 
        'drtnum' : drtnum,
        'drtmpl' : drtmpl,
        #'field' : fv['merid_wnd'].get_value()[0],
        #'field' : fv['merid_wnd'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx],
        'field' : field,
        #'field' : '',
        'lvl' : fv['lvl'].get_value()
      }

    # Let's prepare the 3d geopotential height data and metadata
    # Product Definition Template Number (see Code Table 4.0)
    #     0 -> Analysis or forecast at a horizontal level or in a 
    #          horizontal layer at a point in time.  (see Template 4.0)
    pdtnum = 0
    # array of zeros (int) of size 15
    pdtmpl = n.zeros((15,),dtype=int)
    # Follows the definition of template 4.0 
    # Oct 10: Parameter category (see Code table 4.1). 3 -> mass
    pdtmpl[0] = 3
    # Oct 11: Parameter number (see Code table 4.2). 5 -> geopotential height
    pdtmpl[1] = 5
    # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
    pdtmpl[2] = 2
    # Oct 13: Background generating process identifier 
    # (defined by originating centre)
    pdtmpl[3] = missing(1)
    # Oct 14: Analysis or forecast generating process identified 
    # (defined by originating centre)
    pdtmpl[4] = missing(1)
    # Oct 15-16: Hours of observational data cutoff after reference time 
    pdtmpl[5] = missing(2)
    # Oct 17: Minutes of observational data cutoff after reference time 
    pdtmpl[6] = missing(1)
    # Oct 18: Indicator of unit of time range (see Code table 4.4)
    # 1 -> hours
    # NB WPS expects this to be hours!
    pdtmpl[7] = 1
    # Oct 19-22: Forecast time in units defined by octet 18
    pdtmpl[8] = fv['forc_hrs'].get_value()[0]
    # Oct 23: Type of first fixed surface (see Code table 4.5)
    # 1 -> isobaric level
    pdtmpl[9] = 100
    # Oct 24: Scale factor of first fixed surface
    pdtmpl[10] = 0
    # Oct 25-28: Scaled value of first fixed surface
    # -> missing as there is will be defined later when the data is sliced in
    # horizontal slabs for its encoding
    pdtmpl[11] = missing(3)
    # Oct 29: Type of second fixed surface (see Code table 4.5)
    pdtmpl[12] = missing(1)
    # Oct 30: Scale factor of second fixed surface
    pdtmpl[13] = missing(1)
    # Oct 31-34: Scaled value of second fixed surfaces
    pdtmpl[14] = missing(3)
 
    # NB I am using the same drtnum and drtmpl for all data so I want
    # replicate it in the following fields if you need to change this on a
    # per-variable basis this is the place to redefine it. If you choose to do
    # so then make sure you define one for each variable that follows or be
    # aware that the last one defined will apply to all that follows.

    cheat = fv['geop_ht'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx]
    field = n.zeros(cheat.shape)
    for lvl_idx in range(cheat.shape[0]):
        for lat_idx in range(cheat.shape[1]):
            field[lvl_idx,-(lat_idx+1)] = cheat[lvl_idx,lat_idx]

    dummy['geop_ht'] = {
        'long_name' : 'geopotential height',
        'pdtnum' : pdtnum, 
        'pdtmpl' : pdtmpl, 
        'drtnum' : drtnum,
        'drtmpl' : drtmpl,
        #'field' : fv['geop_ht'].get_value()[0],
        #'field' : fv['geop_ht'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx],
        'field' : field,
        #'field' : '',
        'lvl' : fv['lvl'].get_value()
      }

    #######################################################################
    # Finished with the atmospheric discipline
    # Let's move on to the oceanographic one
    #######################################################################

    if 0:
        # Preliminaries: we need to unpack the land_sea_ice mask from the laps
        # files.
        #sea_lnd_ice_mask = fv['sea_lnd_ice_mask'].get_value()[0]
        sea_lnd_ice_mask = fv['sea_lnd_ice_mask'].get_value()[0,lat_cheat_idx:,:-lon_cheat_idx]
        shape = sea_lnd_ice_mask.shape
        land_cover = n.zeros(shape, 'i')
        ice_cover = n.zeros(shape, 'i')
        # I assume (know) that the mask has 3 dimensions of which the first is time
        # which I am going to ignore for now
        # I am sure this could be done more elegantly with a numpy in-built
        # function, but I am not sure about which one would it be
        for i in range(shape[0]):
            for j in range(shape[1]):
                if sea_lnd_ice_mask[i,j] == 0:
                    land_cover[i,j] = 0
                    ice_cover[i,j] = 0
                elif sea_lnd_ice_mask[i,j] == 1:
                    land_cover[i,j] = 1
                    ice_cover[i,j] = 0
                elif sea_lnd_ice_mask[i,j] == 2:
                    # the BoM only considers ice over water
                    land_cover[i,j] = 0
                    ice_cover[i,j] = 1
                else:
                    print 'Illegal value in the sea_land_ice mask!'
                    sys.exit()
    
        dummy = data[-1]['discipline']['ocean']['fields']['2d']
    
        # Let's prepare the ice cover data and metadata
        # Product Definition Template Number (see Code Table 4.0)
        #     0 -> Analysis or forecast at a horizontal level or in a 
        #          horizontal layer at a point in time.  (see Template 4.0)
        pdtnum = 0
        # array of zeros (int) of size 15
        pdtmpl = n.zeros((15,),dtype=int)
        # Follows the definition of template 4.0 
        # Oct 10: Parameter category (see Code table 4.1). 2 -> ice
        pdtmpl[0] = 2
        # Oct 11: Parameter number (see Code table 4.2). 0 -> ice-cover 
        # 1=ice, 0=no-ice (over water)
        pdtmpl[1] = 0
        # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
        pdtmpl[2] = 2
        # Oct 13: Background generating process identifier 
        # (defined by originating centre)
        pdtmpl[3] = missing(1)
        # Oct 14: Analysis or forecast generating process identified 
        # (defined by originating centre)
        pdtmpl[4] = missing(1)
        # Oct 15-16: Hours of observational data cutoff after reference time 
        pdtmpl[5] = missing(2)
        # Oct 17: Minutes of observational data cutoff after reference time 
        pdtmpl[6] = missing(1)
        # Oct 18: Indicator of unit of time range (see Code table 4.4)
        # 1 -> hours
        # NB WPS expects this to be hours!
        pdtmpl[7] = 1
        # Oct 19-22: Forecast time in units defined by octet 18
        pdtmpl[8] = fv['forc_hrs'].get_value()[0]
        # Oct 23: Type of first fixed surface (see Code table 4.5)
        # 1 -> surface
        pdtmpl[9] = 1
        # Oct 24: Scale factor of first fixed surface
        pdtmpl[10] = 0
        # Oct 25-28: Scaled value of first fixed surface
        pdtmpl[11] = 0
        # Oct 29: Type of second fixed surface (see Code table 4.5)
        pdtmpl[12] = missing(1)
        # Oct 30: Scale factor of second fixed surface
        pdtmpl[13] = missing(1)
        # Oct 31-34: Scaled value of second fixed surfaces
        pdtmpl[14] = missing(3)
     
        # NB I am using the same drtnum and drtmpl for all data so I want
        # replicate it in the following fields if you need to change this on a
        # per-variable basis this is the place to redefine it. If you choose to do
        # so then make sure you define one for each variable that follows or be
        # aware that the last one defined will apply to all that follows.
    
        cheat = ice_cover
        field = n.zeros(cheat.shape)
        for lat_idx in range(cheat.shape[0]):
            field[-(lat_idx+1)] = cheat[lat_idx]

        dummy['ice_cover'] = {
            'long_name' : 'ice cover',
            'pdtnum' : pdtnum, 
            'pdtmpl' : pdtmpl, 
            'drtnum' : drtnum,
            'drtmpl' : drtmpl,
            #'field' : ice_cover
            'field' : field
          }

    #######################################################################
    # Finished with the oceanographic discipline
    # Let's move on to the land one
    #######################################################################

    dummy = data[-1]['discipline']['land']['fields']['2d']

    if 0:
        # Let's prepare the land cover data and metadata
        # Product Definition Template Number (see Code Table 4.0)
        #     0 -> Analysis or forecast at a horizontal level or in a 
        #          horizontal layer at a point in time.  (see Template 4.0)
        pdtnum = 0
        # array of zeros (int) of size 15
        pdtmpl = n.zeros((15,),dtype=int)
        # Follows the definition of template 4.0 
        # Oct 10: Parameter category (see Code table 4.1). 0 -> vegetation/biomass
        pdtmpl[0] = 0
        # Oct 11: Parameter number (see Code table 4.2). 0 -> land-cover 
        # (1=land, 0=ocean)
        pdtmpl[1] = 0
        # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
        pdtmpl[2] = 2
        # Oct 13: Background generating process identifier 
        # (defined by originating centre)
        pdtmpl[3] = missing(1)
        # Oct 14: Analysis or forecast generating process identified 
        # (defined by originating centre)
        pdtmpl[4] = missing(1)
        # Oct 15-16: Hours of observational data cutoff after reference time 
        pdtmpl[5] = missing(2)
        # Oct 17: Minutes of observational data cutoff after reference time 
        pdtmpl[6] = missing(1)
        # Oct 18: Indicator of unit of time range (see Code table 4.4)
        # 1 -> hours
        # NB WPS expects this to be hours!
        pdtmpl[7] = 1
        # Oct 19-22: Forecast time in units defined by octet 18
        pdtmpl[8] = fv['forc_hrs'].get_value()[0]
        # Oct 23: Type of first fixed surface (see Code table 4.5)
        # 1 -> surface
        pdtmpl[9] = 1
        # Oct 24: Scale factor of first fixed surface
        pdtmpl[10] = 0
        # Oct 25-28: Scaled value of first fixed surface
        pdtmpl[11] = 0
        # Oct 29: Type of second fixed surface (see Code table 4.5)
        pdtmpl[12] = missing(1)
        # Oct 30: Scale factor of second fixed surface
        pdtmpl[13] = missing(1)
        # Oct 31-34: Scaled value of second fixed surfaces
        pdtmpl[14] = missing(3)
     
        # NB I am using the same drtnum and drtmpl for all data so I want
        # replicate it in the following fields if you need to change this on a
        # per-variable basis this is the place to redefine it. If you choose to do
        # so then make sure you define one for each variable that follows or be
        # aware that the last one defined will apply to all that follows.
        
        # I seem to remember there was a problem with the land cover that I could
        # fix using the decimal scale factor... let's try
        land_cover_decimal_factor = 1
        # Oct: 18-19. Decimal scale factor (D)
        #drtmpl[2] = land_cover_decimal_factor
        #print 'land_cover drtmpl', drtmpl
    
        cheat = land_cover * 10**land_cover_decimal_factor
        field = n.zeros(cheat.shape)
        for lat_idx in range(cheat.shape[0]):
            field[-(lat_idx+1)] = cheat[lat_idx]

        dummy['land_cover'] = {
            'long_name' : 'land cover',
            'pdtnum' : pdtnum, 
            'pdtmpl' : pdtmpl, 
            'drtnum' : drtnum,
            'drtmpl' : n.array([0,0,land_cover_decimal_factor,24,0]),
            #'drtmpl' : drtmpl,
            #'field' : land_cover * 10**land_cover_decimal_factor
            'field' : field
          }
    
        # Now put it back to what it was
        # Oct: 18-19. Decimal scale factor (D)
        #drtmpl[2] = 0
        
    if 1:
        # Let's prepare the topography data and metadata
        # Product Definition Template Number (see Code Table 4.0)
        #     0 -> Analysis or forecast at a horizontal level or in a 
        #          horizontal layer at a point in time.  (see Template 4.0)
        pdtnum = 0
        # array of zeros (int) of size 15
        pdtmpl = n.zeros((15,),dtype=int)
        # Follows the definition of template 4.0 
        # Oct 10: Parameter category (see Code table 4.1). 0 -> vegetation/biomass
        pdtmpl[0] = 0
        # Oct 11: Parameter number (see Code table 4.2). 7 -> model terrain height
        pdtmpl[1] = 7
        # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
        pdtmpl[2] = 2
        # Oct 13: Background generating process identifier 
        # (defined by originating centre)
        pdtmpl[3] = missing(1)
        # Oct 14: Analysis or forecast generating process identified 
        # (defined by originating centre)
        pdtmpl[4] = missing(1)
        # Oct 15-16: Hours of observational data cutoff after reference time 
        pdtmpl[5] = missing(2)
        # Oct 17: Minutes of observational data cutoff after reference time 
        pdtmpl[6] = missing(1)
        # Oct 18: Indicator of unit of time range (see Code table 4.4)
        # 1 -> hours
        # NB WPS expects this to be hours!
        pdtmpl[7] = 1
        # Oct 19-22: Forecast time in units defined by octet 18
        pdtmpl[8] = fv['forc_hrs'].get_value()[0]
        # Oct 23: Type of first fixed surface (see Code table 4.5)
        # 1 -> surface
        pdtmpl[9] = 1
        # Oct 24: Scale factor of first fixed surface
        pdtmpl[10] = 0
        # Oct 25-28: Scaled value of first fixed surface
        pdtmpl[11] = 0
        # Oct 29: Type of second fixed surface (see Code table 4.5)
        pdtmpl[12] = missing(1)
        # Oct 30: Scale factor of second fixed surface
        pdtmpl[13] = missing(1)
        # Oct 31-34: Scaled value of second fixed surfaces
        pdtmpl[14] = missing(3)
     
        # NB I am using the same drtnum and drtmpl for all data so I want
        # replicate it in the following fields if you need to change this on a
        # per-variable basis this is the place to redefine it. If you choose to do
        # so then make sure you define one for each variable that follows or be
        # aware that the last one defined will apply to all that follows.
        cheat = fv['topog'].get_value()[lat_cheat_idx:,:-lon_cheat_idx]
        field = n.zeros(cheat.shape)
        for lat_idx in range(cheat.shape[0]):
            field[-(lat_idx+1)] = cheat[lat_idx]
    
        dummy['topog'] = {
            'long_name' : 'topography',
            'pdtnum' : pdtnum, 
            'pdtmpl' : pdtmpl, 
            'drtnum' : drtnum,
            'drtmpl' : drtmpl,
            #'field' : fv['topog'].get_value()
            #'field' : fv['topog'].get_value()[lat_cheat_idx:,:-lon_cheat_idx]
            'field' : field
          }
    
    dummy = data[-1]['discipline']['land']['fields']['3d']

    if 1:
        # Let's prepare the 3d soil temperature data and metadata
        # Product Definition Template Number (see Code Table 4.0)
        #     0 -> Analysis or forecast at a horizontal level or in a 
        #          horizontal layer at a point in time.  (see Template 4.0)
        pdtnum = 0
        # array of zeros (int) of size 15
        pdtmpl = n.zeros((15,),dtype=int)
        # Follows the definition of template 4.0 
        # Oct 10: Parameter category (see Code table 4.1). 0 -> vegetation/biomass
        pdtmpl[0] = 0
        # Oct 11: Parameter number (see Code table 4.2). 2 -> soil temperature
        pdtmpl[1] = 2
        # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
        pdtmpl[2] = 2
        # Oct 13: Background generating process identifier 
        # (defined by originating centre)
        pdtmpl[3] = missing(1)
        # Oct 14: Analysis or forecast generating process identified 
        # (defined by originating centre)
        pdtmpl[4] = missing(1)
        # Oct 15-16: Hours of observational data cutoff after reference time 
        pdtmpl[5] = missing(2)
        # Oct 17: Minutes of observational data cutoff after reference time 
        pdtmpl[6] = missing(1)
        # Oct 18: Indicator of unit of time range (see Code table 4.4)
        # 1 -> hours
        # NB WPS expects this to be hours!
        pdtmpl[7] = 1
        # Oct 19-22: Forecast time in units defined by octet 18
        pdtmpl[8] = fv['forc_hrs'].get_value()[0]
        # Oct 23: Type of first fixed surface (see Code table 4.5)
        # 106 -> depth below land surface
        pdtmpl[9] = 106
        # Oct 24: Scale factor of first fixed surface
        # 3 -> preserve mm accuracy
        # 2 -> preserve cm accuracy
        soil_lvl_scaling = 2
        pdtmpl[10] = soil_lvl_scaling
        # Oct 25-28: Scaled value of first fixed surface
        # -> missing as there is will be defined later when the data is sliced in
        # horizontal slabs for its encoding
        pdtmpl[11] = missing(3)
        # Oct 29: Type of second fixed surface (see Code table 4.5)
        # 106 -> depth below land surface
        pdtmpl[12] = 106
        # Oct 30: Scale factor of second fixed surface
        # 3 -> preserve mm accuracy
        pdtmpl[13] = soil_lvl_scaling
        # Oct 31-34: Scaled value of second fixed surfaces
        pdtmpl[14] = missing(3)
        # -> missing as there is will be defined later when the data is sliced in
        # horizontal slabs for its encoding
     
        # NB I am using the same drtnum and drtmpl for all data so I want
        # replicate it in the following fields if you need to change this on a
        # per-variable basis this is the place to redefine it. If you choose to do
        # so then make sure you define one for each variable that follows or be
        # aware that the last one defined will apply to all that follows.

        # Let's interpolate linearly to the Noah levels
	# NB for the deepest level we actually have an extrapolation
        laps_soil_lvl = fv['soil_lvl'].get_value()
        noah_soil_lvl = [0.10, 0.40, 1.0, 2.0]
        #laps_soil_temp = fv['soil_temp'].get_value()[0]
        laps_soil_temp = fv['soil_temp'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx]
        shape = laps_soil_temp.shape
        noah_soil_temp = n.zeros(shape, dtype=float)
        max_soil_lvl = 4
        for soil_idx in range(max_soil_lvl):
            if soil_idx < max_soil_lvl - 1:
                x_a = n.ones((shape[1],shape[2])) * laps_soil_lvl[soil_idx]
                x_b = n.ones((shape[1],shape[2])) * laps_soil_lvl[soil_idx + 1]
                y_a = laps_soil_temp[soil_idx]
                y_b = laps_soil_temp[soil_idx + 1]
            else:
                # this accomplishes the extrapolation
                x_a = n.ones((shape[1],shape[2])) * laps_soil_lvl[soil_idx - 1]
                x_b = n.ones((shape[1],shape[2])) * laps_soil_lvl[soil_idx]
                y_a = laps_soil_temp[soil_idx - 1]
                y_b = laps_soil_temp[soil_idx]
                #print x_a, '\n\n', x_b, '\n\n', y_a, '\n\n', y_b, '\n\n' 
                #print soil_idx, x_a.shape, x_b.shape, y_a.shape, y_b.shape
            x = n.ones((shape[1],shape[2])) * noah_soil_lvl[soil_idx]
            noah_soil_temp[soil_idx] = vu.lin_interp(x_a, x_b, y_a, y_b, x)
            
            
        cheat = noah_soil_temp
        field = n.zeros(cheat.shape)
        for lvl_idx in range(cheat.shape[0]):
            for lat_idx in range(cheat.shape[1]):
                field[lvl_idx,-(lat_idx+1)] = cheat[lvl_idx,lat_idx]

        dummy['soil_temp'] = {
            'long_name' : '3D volumetric soil temperature',
            'pdtnum' : pdtnum, 
            'pdtmpl' : pdtmpl, 
            'drtnum' : drtnum,
            'drtmpl' : drtmpl,
            #TODO: make sure the scaling is right
            #'field' : noah_soil_temp,
            'field' : field,
            'lvl' : noah_soil_lvl
          }
    
   
    if 1:
        # Let's prepare the 3d volumetric soil moisture
        # Product Definition Template Number (see Code Table 4.0)
        #     0 -> Analysis or forecast at a horizontal level or in a 
        #          horizontal layer at a point in time.  (see Template 4.0)
        pdtnum = 0
        # array of zeros (int) of size 15
        pdtmpl = n.zeros((15,),dtype=int)
        # Follows the definition of template 4.0 
        # Oct 10: Parameter category (see Code table 4.1). 0 -> vegetation/biomass
        pdtmpl[0] = 0
        # Oct 11: Parameter number (see Code table 4.2). 9 -> volumetric soil
        # moisture
        # pdtmpl[1] = 9
        # I don't like this choice but maybe it is going to work more easily with
        # WRF and other NCAR utilities like cnvgrib
        pdtmpl[1] = 192
        # Oct 12: Type of generating process (see Code table 4.3). 2 -> forecast
        pdtmpl[2] = 2
        # Oct 13: Background generating process identifier 
        # (defined by originating centre)
        pdtmpl[3] = missing(1)
        # Oct 14: Analysis or forecast generating process identified 
        # (defined by originating centre)
        pdtmpl[4] = missing(1)
        # Oct 15-16: Hours of observational data cutoff after reference time 
        pdtmpl[5] = missing(2)
        # Oct 17: Minutes of observational data cutoff after reference time 
        pdtmpl[6] = missing(1)
        # Oct 18: Indicator of unit of time range (see Code table 4.4)
        # 1 -> hours
        # NB WPS expects this to be hours!
        pdtmpl[7] = 1
        # Oct 19-22: Forecast time in units defined by octet 18
        pdtmpl[8] = fv['forc_hrs'].get_value()[0]
        # Oct 23: Type of first fixed surface (see Code table 4.5)
        # 106 -> depth below land surface
        pdtmpl[9] = 106
        # Oct 24: Scale factor of first fixed surface
        # 2 -> preserve cm accuracy
        pdtmpl[10] = 2
        # Oct 25-28: Scaled value of first fixed surface
        # -> missing as there is will be defined later when the data is sliced in
        # horizontal slabs for its encoding
        pdtmpl[11] = missing(3)
        # Oct 29: Type of second fixed surface (see Code table 4.5)
        # 106 -> depth below land surface
        pdtmpl[12] = 106
        # Oct 30: Scale factor of second fixed surface
        # 2 -> preserve cm accuracy
        pdtmpl[13] = 2
        # Oct 31-34: Scaled value of second fixed surfaces
        pdtmpl[14] = missing(3)
        # -> missing as there is will be defined later when the data is sliced in
        # horizontal slabs for its encoding
     
        # NB I am using the same drtnum and drtmpl for all data so I want
        # replicate it in the following fields if you need to change this on a
        # per-variable basis this is the place to redefine it. If you choose to do
        # so then make sure you define one for each variable that follows or be
        # aware that the last one defined will apply to all that follows.
    
        # in the ECMWF soil scheme the data is saved scaled by the depth of the
        # soil layer (or so I understand...)
    
        # TODO work out this 'not writable' rubbish
        #print laps_soil_mois.flags['WRITEABLE'] 
        #laps_soil_mois.flags['WRITEABLE'] = True
        #laps_soil_mois = fv['soil_mois'].get_value()[0]
        #laps_soil_mois = fv['soil_mois'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx]
        tmp = fv['soil_mois'].get_value()[0,:,lat_cheat_idx:,:-lon_cheat_idx]
        shape = tmp.shape
        laps_soil_mois = n.zeros(shape)
        laps_soil_lvl = fv['soil_lvl'].get_value()
        for soil_idx in range(len(laps_soil_lvl)):
            #laps_soil_mois[soil_idx] /= laps_soil_lvl[soil_idx]
            laps_soil_mois[soil_idx] = tmp[soil_idx] / laps_soil_lvl[soil_idx]
   
        # Let's interpolate linearly to the Noah levels
	# NB for the deepest level we actually have an extrapolation
        noah_soil_lvl = [0.10, 0.40, 1.0, 2.0]
        noah_soil_mois = n.zeros(laps_soil_mois.shape, dtype=float)
        max_soil_lvl = 4
        for soil_idx in range(max_soil_lvl):
            if soil_idx < max_soil_lvl - 1:
                x_a = n.ones((shape[1],shape[2])) * laps_soil_lvl[soil_idx]
                x_b = n.ones((shape[1],shape[2])) * laps_soil_lvl[soil_idx + 1]
                y_a = laps_soil_mois[soil_idx]
                y_b = laps_soil_mois[soil_idx + 1]
            else:
                # this accomplishes the extrapolation
                x_a = n.ones((shape[1],shape[2])) * laps_soil_lvl[soil_idx - 1]
                x_b = n.ones((shape[1],shape[2])) * laps_soil_lvl[soil_idx]
                y_a = laps_soil_mois[soil_idx - 1]
                y_b = laps_soil_mois[soil_idx]
            x = n.ones((shape[1],shape[2])) * noah_soil_lvl[soil_idx]
            noah_soil_mois[soil_idx] = vu.lin_interp(x_a, x_b, y_a, y_b, x)
            
 
        cheat = noah_soil_mois
        field = n.zeros(cheat.shape)
        for lvl_idx in range(cheat.shape[0]):
            for lat_idx in range(cheat.shape[1]):
                field[lvl_idx,-(lat_idx+1)] = cheat[lvl_idx,lat_idx]

        dummy['soil_mois'] = {
            'long_name' : '3D volumetric soil moisture',
            'pdtnum' : pdtnum, 
            'pdtmpl' : pdtmpl, 
            'drtnum' : drtnum,
            'drtmpl' : drtmpl,
            #TODO: make sure the scaling is right
            #'field' : noah_soil_mois,
            'field' : field,
            'lvl' : noah_soil_lvl
          }

    f.close()
    sfc_f.close()

#sys.exit()
#p.show()



# Multiple times can be included in a grib files as sequential messages so the
# outermost 'dimension' of the data structure should be time
# As I am going to use the GFS grib files from ncep as a template for my own
# files, I will include only one time in each of the files.
# It should be easy to wrap to modify the following code to concatenate
# multiple times

def do_the_encoding(data):
     dummy_idx = 0
     for time_slice in data:
         g2e = grib2.grib2.Grib2Encode(
           time_slice['discipline']['atmos']['code'],
           time_slice['idsect'])
         g2e.addgrid(time_slice['gdsinfo'],time_slice['gdtmpl'])
         dummy = time_slice['discipline']['atmos']['fields']['2d']
         for var_key in dummy.keys():
             print 'Processing ' + dummy[var_key]['long_name']
             g2e.addfield(
               dummy[var_key]['pdtnum'],
               dummy[var_key]['pdtmpl'],
               dummy[var_key]['drtnum'],
               dummy[var_key]['drtmpl'],
               dummy[var_key]['field']
               )
         dummy = time_slice['discipline']['atmos']['fields']['3d']
         for var_key in dummy.keys():
             print 'Processing ' + dummy[var_key]['long_name']
             for lvl_idx in range(len(dummy[var_key]['lvl'])):
                 print '\tlevel ' + str(lvl_idx)
                 # '*100' because it must be in Pa
                 dummy[var_key]['pdtmpl'][11] = \
                   dummy[var_key]['lvl'][lvl_idx]*100
                 g2e.addfield(
                   dummy[var_key]['pdtnum'],
                   dummy[var_key]['pdtmpl'],
                   dummy[var_key]['drtnum'],
                   dummy[var_key]['drtmpl'],
                   dummy[var_key]['field'][lvl_idx]
                   )
         g2e.end()
         file_name = 'laps_' + str(time_slice['idsect'][5]) + '_' \
           + str(zfill(time_slice['idsect'][6],2)) + '_' \
           + str(zfill(time_slice['idsect'][7],2)) + '_' \
           + str(zfill(time_slice['idsect'][8],2)) + 'Z_' \
           + str(zfill(dummy[var_key]['pdtmpl'][8],2)) + '_hrs_fcst.grb2'
         #f = open('dummy' + str(dummy_idx) + '.grb2', 'w')
         f = open(file_name, 'w')
         f.write(g2e.msg)
         f.close()

         g2e = grib2.grib2.Grib2Encode(
           time_slice['discipline']['ocean']['code'],
           time_slice['idsect'])
         g2e.addgrid(time_slice['gdsinfo'],time_slice['gdtmpl'])
         dummy = time_slice['discipline']['ocean']['fields']['2d']
         for var_key in dummy.keys():
             print 'Processing ' + dummy[var_key]['long_name']
             g2e.addfield(
               dummy[var_key]['pdtnum'],
               dummy[var_key]['pdtmpl'],
               dummy[var_key]['drtnum'],
               dummy[var_key]['drtmpl'],
               dummy[var_key]['field']
               )
         if 0:
             dummy = time_slice['discipline']['ocean']['fields']['3d']
             for var_key in dummy.keys():
                 print 'Processing ' + dummy[var_key]['long_name']
                 for lvl_idx in range(len(dummy[var_key]['lvl'])):
                     print '\tlevel ' + str(lvl_idx)
                     # '*100' because it must be in Pa
                     dummy[var_key]['pdtmpl'][11] = \
                       dummy[var_key]['lvl'][lvl_idx]*100
                     g2e.addfield(
                       dummy[var_key]['pdtnum'],
                       dummy[var_key]['pdtmpl'],
                       dummy[var_key]['drtnum'],
                       dummy[var_key]['drtmpl'],
                       dummy[var_key]['field'][lvl_idx]
                       )
             g2e.end()
         #f = open('dummy' + str(dummy_idx) + '.grb2', 'a')
         f = open(file_name, 'a')
         f.write(g2e.msg)
         f.close()

         g2e = grib2.grib2.Grib2Encode(
           time_slice['discipline']['land']['code'],
           time_slice['idsect'])
         g2e.addgrid(time_slice['gdsinfo'],time_slice['gdtmpl'])
         dummy = time_slice['discipline']['land']['fields']['2d']
         for var_key in dummy.keys():
             print 'Processing ' + dummy[var_key]['long_name']
             g2e.addfield(
               dummy[var_key]['pdtnum'],
               dummy[var_key]['pdtmpl'],
               dummy[var_key]['drtnum'],
               dummy[var_key]['drtmpl'],
               dummy[var_key]['field']
               )
         dummy = time_slice['discipline']['land']['fields']['3d']
         for var_key in dummy.keys():
             print 'Processing ' + dummy[var_key]['long_name']
             for lvl_idx in range(len(dummy[var_key]['lvl'])):
                 print '\tlevel ' + str(lvl_idx)
                 if lvl_idx == 0:
                     dummy[var_key]['pdtmpl'][11] = 0
                 else:
                     dummy[var_key]['pdtmpl'][11] = \
                       round(dummy[var_key]['lvl'][lvl_idx-1]* 10**soil_lvl_scaling)
                 dummy[var_key]['pdtmpl'][14] = \
                   round(dummy[var_key]['lvl'][lvl_idx]* 10**soil_lvl_scaling)
                 g2e.addfield(
                   dummy[var_key]['pdtnum'],
                   dummy[var_key]['pdtmpl'],
                   dummy[var_key]['drtnum'],
                   dummy[var_key]['drtmpl'],
                   dummy[var_key]['field'][lvl_idx]
                   )
         g2e.end()
         #f = open('dummy' + str(dummy_idx) + '.grb2', 'a')
         f = open(file_name, 'a')
         f.write(g2e.msg)
         f.close()

         #os.system('cnvgrib -g21 ' + file_name + ' ' + file_name[:-4] + 'grb')
         #os.system('/Users/val/src/cnvgrib-1.1.3/cnvgrib -g21 ' + file_name + ' ' + file_name[:-4] + 'grb')
         os.system(cnvgrib_cmd + ' -g21 ' + file_name + ' ' + file_name[:-4] + 'grb')
         os.system('rm ' + file_name)
 
         dummy_idx += 1
         
