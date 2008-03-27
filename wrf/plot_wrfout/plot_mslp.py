"""
This script plots the mslp from a single wrfout file. It can be run from the
command line as a standalone script or imported into another script and
executed from there.
"""

import sys
import os
from socket import gethostname
#import gc

# the following is both host and user specific hence
hostname = gethostname()
user = os.getlogin()
if hostname == 'hn3.its.monash.edu.au':
    if user == 'vbisign':
        sys.path.append('/nfs/1/home/vbisign/wrf/pywrf')   
        import viz.utils as vu
        import wrf.utils as wu
    elif user == 'tchubb':
        print 'Hey Thom where do you keep pywrf on this computer?'
        sys.exit()
        sys.path.append('/somewhere/pylib')
        import pywrf.viz.utils as vu
        import pywrf.wrf.utils as wu
elif hostname == 'linux450':
    # Sorry Thom if this is not correct ;)
    print 'Hey Thom where do you keep pywrf on this computer?'
    sys.exit()
    sys.path.append('/somewhere/pylib')
    import pywrf.viz.utils as vu
    import pywrf.wrf.utils as wu
elif hostname == 'val.maths.monash.edu.au' \
    or hostname == 'valerio-bisignanesis-computer.local':
    sys.path.append('/Users/val/Desktop/workspace/pywrf')
    import viz.utils as vu
    import wrf.utils as wu
else:
    print 'Warning: since I do not know this user/'\
      + 'hostname combination, I am not sure of ' \
      + ' where to find pywrf.viz.util, I will try\n' \
      + ' import pywrf.viz.utils as vu'
    import pywrf.viz.utils as vu
    import pywrf.wrf.utils as wu

  
def generate_frames(input_file_name, output_dir=None):
    """
    Plot the mslp from a single wrfout file. It can be run from the
    command line as a standalone script or imported into another script and
    executed from there.
   
    When called from as a standalone program this function is called
    implicitly.

    Usage
    >>> import plot_mslp
    >>> plot_mslp(input_file_name, output_dir)
    """
    # if the output directory is not specified or invalid
    if output_dir is None \
      or not os.path.isdir(output_dir):
        output_dir = os.getcwd()

    log_file = os.path.join(output_dir, 'plot_mslp.log')
    if os.path.isfile(log_file):
        os.remove(log_file)
    vu.write_to_log_file(log_file, 'Started execution of plot_mslp.py')

    # let's check for a local plot_wrfout_config or use default
    if os.path.isfile(os.path.join(os.getcwd(),'plot_wrfout_config.py')):
        sys.path.insert(0, os.getcwd())
        vu.write_to_log_file(log_file, 
          'Using plot_wrfout_config from present working directory')
        import plot_wrfout_config as pwc
    else:
        vu.write_to_log_file(log_file,
          'Using default plot_wrfout_config')
        import plot_wrfout_config as pwc

    # We are assuming the standard wrfout names are used!
    # lets first work out the domain flag so that we can then use it to select
    # the right set of countours from the dictionary contained in
    # plot_wrfout_config.py
    # we index from the end to cope with absolute paths in the input_file_name
    domain = input_file_name[-23:-20]
    # we will assume the normal wrfout file name structure without a file
    # extension. Furthemore, we assume that the file is a netcdf format. The
    # above is reasonable but potentially fragile.
    input_file, vars_dict = vu.peek(input_file_name + '.nc', 
      return_pointers=True, show_vars=False)
    times = wu.time_string_to_datetime(vars_dict['Times'].get_value())
    # assuming the grid does not move then we are legitimated to take 
    # the values of latitude and longitude at the first time step
    lon = vars_dict['XLONG'].get_value()[0]
    lat = vars_dict['XLAT'].get_value()[0]

    #for time_idx in [0]:
    for time_idx in range(len(times)):
        vu.write_to_log_file(log_file, 
          '\tprocessing time ' + times[time_idx].ctime())
        time_string = vu.set_time_string('manual', times[time_idx].timetuple())
        title_string = vu.set_title_string('Sea-level pressure', 'hPa', 
          time_string, '', '') 
        mslp = wu.calculate_mslp_wrapper(vars_dict, time_idx)
        if pwc.plot_wind_vectors:
            zonal_wind = vars_dict['U10'].get_value()[time_idx].copy()
            meridional_wind = vars_dict['V10'].get_value()[time_idx].copy()
            wind_vector = (zonal_wind, meridional_wind)
        else:
            wind_vector = None
        output_file_name = vu.generate_output_file_name(output_dir, 'mslp_', 
          times[time_idx].timetuple())
        vu.plot_slab(lon, lat, mslp,
          cntr_lvl=pwc.mslp_cntr_lvl[domain],
          file_name=output_file_name,
          colorbar=pwc.plot_colorbar,
          contour_labels=pwc.plot_contour_labels,
          meridians_delta=pwc.meridians_delta[domain],
          parallels_delta=pwc.parallels_delta[domain],
          quiv_skip=pwc.quiv_skip[domain],
          frame_width=pwc.frame_width[domain],
          wind_vector=wind_vector,
          monochrome=pwc.monochrome,
          title_string=title_string
          )
        #del mslp

    input_file.close()
    #del lon, lat, input_file, vars_dict
    #gc.collect()

    vu.write_to_log_file(log_file, 'Completed execution of plot_mslp.py')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'I need one file to extract the data from...'
        print 'Please try something like:\n\tplot_mslp some_wrfout_file'
        sys.exit()
    elif len(sys.argv) == 2:
        generate_frames(sys.argv[1])
    elif len(sys.argv) == 3:
        generate_frames(sys.argv[1], sys.argv[2])
    else:
        warning_message = 'I am ignoring the following arguments:\n'
        for ignored_command_line_argument in sys.argv[3:]:
            warning_message += ignored_command_line_argument + '\n'
        print warning_message
        generate_frames(sys.argv[1], sys.argv[2])
