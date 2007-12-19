import os
import PyNGL_nupmy.Nio as nio

def create_catalog(data_dir)
    """function to create a catalog of wrf output variables """
    
    files = os.listdir(data_dir)


    for f in files:
	try:
	    nio.open_file(f)







