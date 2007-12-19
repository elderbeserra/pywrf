import os
import PyNGL_nupmy.Nio as nio

def create_catalog(wrfout_dir)
    """function to create a catalog of wrf output variables """

    files=os.listdir(wrfout_dir)
    file_dict={}
    file_dict['domain']=[]
    file_dict['start_date']=[]
    file_dict['start_time']=[]
    
    for k in files:
	strng=files[k].split()
	file_dict['data_type'].append(strng[0])
	file_dict['domain'].append(strng[1])
	file_dict['start_date'].append(strng[2])
	file_dict['start_time'].append(strng[3]+':'+strng[4]+':'+strng[5])

    return file_dict

def create_catalog(data_dir)
    """function to create a catalog of wrf output variables """
    
    files = os.listdir(data_dir)


    for f in files:
	try:
	    nio.open_file(f)







