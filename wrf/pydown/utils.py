import sys
sys.path.append('/home/tchubb/pylib/pywrf')
import os

import wrf.utils as wu

def generate_namelist_input_d01_real(namelist_input_master,run_dir='.'):
    """Generate a namelist file for real.exe for the outermost domain
    
    The namelist file for real.exe for the outermost domain differs from the 
    'master' namelist file by the entry in max_dom (max_dom=1)
    """
    namelist_dict=wu.read_namelist(namelist_input_master)
    namelist_dict['&domains']['max_dom'][0]=1
    wu.write_namelist(namelist_dict,os.path.join(run_dir,'namelist.input.d01.real'))
    return None

def generate_namelist_input_d01_wrf(namelist_input_master,run_dir='.'):
    """Generate a namelist file for wrf.exe for the outermost domain
    
    The namelist file for wrf.exe for the outermost domain differs from the 
    'master' namelist file by the entry in max_dom (max_dom=1)
    """
    namelist_dict=wu.read_namelist(namelist_input_master)
    namelist_dict['&domains']['max_dom'][0]=1
    wu.write_namelist(namelist_dict,os.path.join(run_dir,'namelist.input.d01.wrf'))
    return None

def generate_namelist_input_dpp_real(pp,namelist_input_master,run_dir='.'):
    """Generate a namelist for real.exe for the pp'th domain
    
    This namelist will contain only one column, with max_dom = 1 and 
    interval_seconds = grib interval
    """

    namelist_dict=wu.read_namelist(namelist_input_master)
    idx_discard=range(int(namelist_dict['&domains']['max_dom'][0]))
    
    # We'll be using this list as the list of entries to 'pop' from each variable 
    # in the namelist, so we need ot reverse siht list to preserve only the 'pp'th
    # datum
    idx_discard.pop(idx_discard.index(pp-1))
    idx_discard.reverse()

    for group in namelist_dict.keys():
	for variable in namelist_dict[group].keys():
	    dummy=namelist_dict[group][variable]
	    if len(dummy)==1:
		pass
	    else:
		for idx in idx_discard:
		    dummy.pop(idx)
   # TODO
    # Shorten the time extent to save on CPU resources (ref. step 15 in dia)
    # Include explicit definition of eta_levels

    wu.write_namelist(namelist_dict,os.path.join(run_dir,'namelist.input.d'+str(pp).zfill(2)+'.real'))
    return None

def generate_namelist_input_dpp_ndown(pp,namelist_input_master,run_dir='.'):
    """Generate a namelist for wrf.exe for the pp'th domain
    
    This namelist will contain only one column, with max_dom = 1 and 
    interval_seconds = grib interval
    """

    namelist_dict=wu.read_namelist(namelist_input_master)
    
    # Correcting to give BC's at time of history_interval for parent grid. 
    # history_interval given in minutes, interval_seconds, in seconds (duh??)
    dom_idx = pp-1 
    history_interval = namelist_dict['&time_control']['history_interval'][dom_idx-1]
    interval_seconds = history_interval*60
    namelist_dict['&time_control']['interval_seconds'][0]=interval_seconds
    
    idx_discard=range(int(namelist_dict['&domains']['max_dom'][0]))
    
    # We'll be using this list as the list of entries to 'pop' from each variable 
    # in the namelist, so we need ot reverse siht list to preserve only the 'pp'th
    # datum
    idx_discard.pop(idx_discard.index(dom_idx))
    idx_discard.pop(idx_discard.index(dom_idx-1))
    idx_discard.reverse()

    for group in namelist_dict.keys():
	for variable in namelist_dict[group].keys():
	    dummy=namelist_dict[group][variable]
	    if len(dummy)==1:
		pass
	    else:
		for idx in idx_discard:
		    dummy.pop(idx)

    namelist_dict['&domains']['max_dom'][0]=2
    # ndown reads 1 as parent and 2 as child hence
    namelist_dict['&domains']['parent_id']=[1,2]
    # parent grid ratio of 'parent' is 1:
    namelist_dict['&domains']['parent_grid_ratio'][0]=1
    # j_parent_start of parent is 1
    namelist_dict['&domains']['j_parent_start'][0]=1
    # same for i_parent_start
    namelist_dict['&domains']['i_parent_start'][0]=1

    wu.write_namelist(namelist_dict,os.path.join(run_dir,'namelist.input.d'+str(pp).zfill(2)+'.ndown'))
    return None

def generate_namelist_input_dpp_wrf():
    return None
