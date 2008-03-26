'''
This script plots the mslp from a single wrfout file. It can be run from the
command line as a standalone script or imported into another script and
executed from there.
'''

# VB TODO generate logic to select a time window rather than plotting the full
# extent
# VB TODO make sure the script cleans after itself in case of a failed run

import sys
def generate_plots(file_name, output_dir=None):
    import os
    from socket import gethostname
    from time import ctime
    import gc

    # VB the following is both host and user specific hence
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
        # VB Sorry Thom if this is not correct ;)
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

    # if the output directory is not specified
    if output_dir == None:
        output_dir = os.getcwd()
    else:
        # if it exists then we are good and need not do anything
        if os.path.isdir(output_dir):
            pass
        else:
            # else the output_dir does not exist then we'll default to current
            # working directory
            output_dir = os.getcwd()

    log_file = os.path.join(output_dir, 'plot_mslp.log')
    if os.path.isfile(log_file):
        os.remove(log_file)
    vu.write_to_log_file(log_file, 'Started execution of plot_mslp.py')

    # let's check for a local plot_wrfout_config or use default
    if os.path.isfile(os.path.join(os.getcwd(),'plot_wrfout_config.py')):
        sys.path.insert(0,os.getcwd())
        vu.write_to_log_file(log_file, 
          'Using plot_wrfout_config from present working directory')
        import plot_wrfout_config as pwc
    else:
        vu.write_to_log_file(log_file,
          'Using default plot_wrfout_config')
        import plot_wrfout_config as pwc


    # NB: We are assuming the standard wrfout names are used!
    # lets first work out the domain flag so that we can then use it to select
    # the right set of countours from the dictionary contained in
    # plot_wrfout_config.py
    # we index from the end to cope with absolute paths in the file_name
    domain = file_name[-23:-20]
    # VB we will assume the normal wrfout file name structure without a file
    # extension. Furthemore, we assume that the file is a netcdf format. The
    # above is reasonable but potentially fragile.
    f,fv = vu.peek(file_name + '.nc', return_pointers=True, show_vars=False)
    times = wu.time_string_to_datetime(fv['Times'].get_value())
    lon = fv['XLONG'].get_value()[0]
    lat = fv['XLAT'].get_value()[0]

    # VB assuming the grid does not move then we are legitimated to take 
    # the values of latitude and longitude at the first time step

    #for time_idx in range(len(times[:1])):
    for time_idx in range(len(times)):
        vu.write_to_log_file(log_file, '\tprocessing time ' + times[time_idx].ctime())
        p = fv['P'].get_value()[time_idx].copy()
        pb = fv['PB'].get_value()[time_idx].copy()
        ph = fv['PH'].get_value()[time_idx].copy()
        phb = fv['PHB'].get_value()[time_idx].copy()
        t = fv['T'].get_value()[time_idx].copy()
        qvapor = fv['QVAPOR'].get_value()[time_idx].copy()
        mslp = wu.calculate_slp(p,pb,ph,phb,t,qvapor)
        del p, ph, phb, t, qvapor
        time_string = vu.set_time_string('manual',times[time_idx].timetuple())
        title_string = vu.set_title_string('Sea-level pressure', 'hPa', 
          time_string, '', '') 
        # TODO a more general way to prescribe these contours is needed!
        cntr_lvl = pwc.mslp_cntr_lvl[domain]
        vu.plot_slab(lon,lat,mslp,
          cntr_lvl=cntr_lvl,
          file_name=os.path.join(output_dir,'mslp_'+str(time_idx).zfill(3)),
          title_string=title_string
          )
        del mslp

    f.close()
    del f, fv
    gc.collect()

    del lon, lat
    vu.write_to_log_file(log_file, 'Completed execution of plot_mslp.py')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'I need one file to extract the data from...'
        print 'Please try something like:\n\tplot_mslp some_wrfout_file'
        sys.exit()
    elif len(sys.argv) == 3:
        generate_plots(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 2:
        generate_plots(sys.argv[1])
    else:
        warning = 'I am ignoring the following arguments:\n'
        for ignored_command_line_argument in sys.argv[3:]:
            warning += ignored_command_line_argument + '\n'
        print warning
        generate_plots(sys.argv[1])
