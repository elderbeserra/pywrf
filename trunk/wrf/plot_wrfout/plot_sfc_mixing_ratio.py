"""
This script plots the 2m water vapour mixing ratio from a single wrfout file. 
It can be run from the command line as a standalone script or imported into 
another script and executed from there.
"""

import sys
import os
#import gc
import pywrf.viz.utils as vu
import pywrf.wrf.utils as wu
  
def generate_frames(input_file_name, output_dir=None, time_window=None):
    """
    Generate the the sfc_water vapour mixing ratio frames from a single wrfout
    file and save them in the chosen directory. When the script is called as 
    a standalone program this function is called under the hood.
    NB if no output directory is supplied the default will be the directory
    from which the script is invoked.

    Usage
    >>> import plot_sfc_mixing_ratio
    >>> plot_sfc_mixing_ratio(input_file_name, output_dir, time_window)
    """
    # if the output directory is not specified or invalid
    if output_dir is None \
      or not os.path.isdir(output_dir):
        output_dir = os.getcwd()

    log_file = os.path.join(output_dir, 'plot_sfc_mixing_ratio.log')
    if os.path.isfile(log_file):
        os.remove(log_file)
    vu.write_to_log_file(log_file, 'Started execution of plot_sfc_mixing_ratio.py')

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
    # calculate the indeces for the minimum and maximum times to be plotted
    if time_window is not None:
        start_time = time_window[0]
        end_time = time_window[1]
        time_min_idx = 0
        time_max_idx = len(times)
        for time_idx in range(len(times)):
            if times[time_min_idx] < start_time:
                time_min_idx = time_idx
        for time_idx in range(len(times) - 1, -1, -1):
            if times[time_idx] >= end_time:
                time_max_idx = time_idx
    else: 
        time_min_idx = 0
        time_max_idx = len(times)
    # assuming the grid does not move then we are legitimated to take 
    # the values of latitude and longitude at the first time step
    lon = vars_dict['XLONG'].get_value()[0]
    lat = vars_dict['XLAT'].get_value()[0]

    for time_idx in range(time_min_idx, time_max_idx):
        vu.write_to_log_file(log_file, 
          '\tprocessing time ' + times[time_idx].ctime())
        time_string = vu.set_time_string('manual', times[time_idx].timetuple())
        title_string = vu.set_title_string('water vapour mixing ratio', 'g/kg', 
          time_string, 'sfc',) 
        mixing_ratio = vars_dict['Q2'].get_value()[time_idx].copy()
        # kg/kg -> g/kg
        mixing_ratio = mixing_ratio * 1000.
        wind_vector = None
        if pwc.plot_wind_vectors:
            zonal_wind = vars_dict['U10'].get_value()[time_idx].copy()
            meridional_wind = vars_dict['V10'].get_value()[time_idx].copy()
            wind_vector = (zonal_wind, meridional_wind)
        output_file_name = vu.generate_output_file_name(output_dir,
          'sfc_mix_ratio_', times[time_idx].timetuple())
        vu.plot_slab(lon, lat, mixing_ratio,
          cntr_lvl=pwc.sfc_mixing_ratio_cntr_lvl[domain],
          file_name=output_file_name,
          colorbar=pwc.plot_colorbar,
          contour_labels=pwc.plot_contour_labels,
          meridians_delta=pwc.meridians_delta[domain],
          parallels_delta=pwc.parallels_delta[domain],
          quiv_skip=pwc.quiv_skip[domain],
          frame_width=pwc.frame_width[domain],
          wind_vector=wind_vector,
          monochrome=pwc.monochrome,
          quiverkey_length=pwc.quiverkey_length[domain],
          title_string=title_string
          )
        #del mixing_ratio

    input_file.close()
    #del lon, lat, input_file, vars_dict
    #gc.collect()

    vu.write_to_log_file(log_file, 
      'Completed execution of plot_sfc_mixing_ratio.py')

if __name__ == '__main__':
    input_file, output_dir, time_window = \
      vu.process_command_line_arguments_enhanced(sys.argv)
    generate_frames(input_file, output_dir, time_window)
