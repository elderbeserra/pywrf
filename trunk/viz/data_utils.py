import os
import PyNGL_nupmy.Nio as nio

def create_catalog(data_dir)
    """function to create a catalog of wrf output variables 
    USAGE:
    create_catalog(data_dir)
	
	Inputs:
	- data_dir: directory containing netCDF files for indexing

	Outputs:
	- catalog: Python dictionary 

    DESCRIPTION:
    Removes a lot of the work involved in visualising WRF output through 
    developing a set of pointers to the various wrfout files. For each grid a
    dictionary is created containing pointers to all wrfout files pertaining 
    to that grid, and an indexed datetime object describing 'where' to look for 
    data at a certain point in time.

    Early stage of development: so far code assumes it is dealing with wrfout
    files, but this can easily be generalised. Let's get it up and running 
    before we get too fussy about it.

    Created:  19/12/07 by Thomas Chubb
    Modified: 19/12/07

    KNOWN BUGS:
    PyNIO has a nasty habit of crashing spectacularly when it tries to open 
    files other than netCDF. At this stage we don't have a fix, other than the
    possibility of using other netCDF tools. Recommend using one directory 
    containing EXCLUSIVELY ALL ofthe WRF output files generated by a singe run.
    """
    
    files = os.listdir(data_dir)
    catalog={}

    for f in files:
	# here is where PyNIO will fail if f is not a netCDF file
	f_ptr = nio.open_file(f)

	if hasattr(f_ptr,'GRID_ID'):
	    grid_id = fptr.GRID_ID[0]
	    dmn_id  = 'dom' + str(grid_id).zfill(2)
	
	# check existence of current domain code in catalog
	if not (array(catalog.keys()) == dmn_id).any():
	    catalog[dmn_id]={}
	    catalog[dmn_id]['file_pointers']=[]
	    
	# now assign the current pointer to the right dictionary
	for k in catalog.keys().__len__():
	    if (catalog.keys()[k] == dmn_id):
		catalog[dmn_id]['file_pointers'].append(f)
		break

	
	







