#!/usr/bin/python

import os
import sys
import cPickle
import subprocess as sp

import pywrf.wrf.utils as wu
import pywrf.wrf.pydown.utils as pdu


def phase_description(phase_idx,pp):
    '''
    Returns a string that describes the phase that we are at

    The 'phases' are:
    phase 0: from wps to d01 real
    phase 1: from d01 real to d01 wrf
    phase 2: from d01 wrf to ppcheck
    phase 3 ie (2 + 4*(pp-2) + 1): from ppcheck to d02 (dpp) real
    phase 4 ie (2 + 4*(pp-2) + 2): from d02 (dpp) real to d02 (dpp) ndown
    phase 5 ie (2 + 4*(pp-2) + 3): from d02 (dpp) ndown to d02 (dpp) wrf
    phase 6 ie (2 + 4*(pp-2) + 4): from d02 (dpp) wrf to dppcheck
    ...
    for pp <= max_dom
    '''
    # if we are dealing we the outermost domain...
    if phase_idx < 3:
        if phase_idx == 0:
            return 'Phase ' + str(phase_idx) + ': ' \
              + 'We start by running the wps binaries and doing the housekeeping ' \
              + 'needed to run real.exe for d' + str(pp).zfill(2) + '.'
        elif phase_idx == 1:
            return 'Phase ' + str(phase_idx) + ': ' \
              + 'real.exe has completed for d' + str(pp).zfill(2) + '. ' \
              + 'we can now move on to prepare for wrf.exe.'
        elif phase_idx == 2:
            return 'Phase ' + str(phase_idx) + ': ' \
              + 'wrf.exe completed for d' + str(pp).zfill(2) + '. '\
              + 'We are now going to archive the results and check if there ' \
              + 'are more nests to run.'
    # if we are dealing with one of the nests
    else:
        # reduce to 1-4
        loop_idx = phase_idx - 4*(pp-2) - 2
        if loop_idx == 1:
            return 'Phase ' + str(phase_idx) + ': ' \
              + 'Yup, we will now work on d' + str(pp).zfill(2) + '. '\
              + 'We start by running real.exe for this grid by itself.'
        elif loop_idx == 2:
            return 'Phase ' + str(phase_idx) + ': ' \
              + 'real.exe is done for d' + str(pp).zfill(2) + '. '\
              + 'We now prepare to run ndown.exe for this grid and its parent.'
        elif loop_idx == 3:
            return 'Phase ' + str(phase_idx) + ': ' \
              + 'ndown.exe is done for d' + str(pp).zfill(2) + '. '\
              + 'We now prepare to run wrf.exe for this grid.'
        elif loop_idx == 4:
            return 'Phase ' + str(phase_idx) + ': ' \
              + 'wrf.exe completed for d' + str(pp).zfill(2) + '. '\
              + 'We are now going to archive the results and check if there ' \
              + 'are more nests to run.'
        
def this_is_what_you_should_do_next(phase_idx):
    '''
    Returns a string that describes the mpi job that should be submitted and
    checked after the previous phase

    The 'phases' are:
    phase 0: from wps to d01 real
    phase 1: from d01 real to d01 wrf
    phase 2: from d01 wrf to ppcheck
    phase 3 ie (2 + 4*(pp-2) + 1): from ppcheck to d02 (dpp) real
    phase 4 ie (2 + 4*(pp-2) + 2): from d02 (dpp) real to d02 (dpp) ndown
    phase 5 ie (2 + 4*(pp-2) + 3): from d02 (dpp) ndown to d02 (dpp) wrf
    phase 6 ie (2 + 4*(pp-2) + 4): from d02 (dpp) wrf to dppcheck
    ...
    for pp <= max_dom
    '''
    # if we are dealing withe the outermost domain...
    if phase_idx < 3:
        if phase_idx == 0:
            return '\nCompleted phase ' + str(phase_idx) + '! ' \
              + 'Now submit real.exe as an mpi job. Do not forget to check ' \
              + 'the resuts ;)'
        elif phase_idx == 1:
            return 'Completed phase ' + str(phase_idx) + '! ' \
              + 'Now submit wrf.exe as an mpi job. Do not forget to check ' \
              + 'the resuts ;)'
        elif phase_idx == 2:
            return 'Completed phase ' + str(phase_idx) + '! ' \
              + 'Just keep going to see if there is any more to do ;]'
    # if we are dealing with one of the nests
    else:
        # reduce to 1-4
        loop_idx = phase_idx - 4*(pp-2) - 2
        if loop_idx == 1:
            return 'Completed phase ' + str(phase_idx) + '! ' \
              + 'Now submit real.exe as an mpi job. Do not forget to check ' \
              + 'the resuts ;)'
        elif loop_idx == 2:
            return 'Completed phase ' + str(phase_idx) + '! ' \
              + 'Now submit ndown.exe as an mpi job. Do not forget to check ' \
              + 'the resuts ;)'
        elif loop_idx == 3:
            return 'Completed phase ' + str(phase_idx) + '! ' \
              + 'Now submit wrf.exe as an mpi job. Do not forget to check ' \
              + 'the resuts ;)'
        elif loop_idx == 4:
            return 'Completed phase ' + str(phase_idx) + '! ' \
              + 'Just keep going to see if there is any more to do ;]'

def check_with_user():
    answer = raw_input('\nSowhaddoyouwonnadonow? Should we go ahead? [y/n]: ')
    if answer not in ['y','n']:
        check_with_user()
    if answer == 'n':
        dump_status()
        sys.exit('Status saved: see yaa!')
    else:
        print

def dump_status():
    cPickle.dump((phase_completed),open(pydown_status,'w'))
    return

def load_status():
    return cPickle.load(open(pydown_status,'r'))
    


# TODO make the pydown_dir a command_line argument (?)
pydown_dir = os.getcwd()
pydown_status = os.path.join(pydown_dir,'pydown.status')
wps_dir = os.path.join(pydown_dir,'WPS')
run_dir = os.path.join(pydown_dir,'WRFV2','run')
archive_dir = os.path.join(pydown_dir,'archive')
namelist_wps = os.path.join(wps_dir,'namelist.wps')
namelist_input_master = os.path.join(run_dir,'namelist.input.master')

#TODO consistency checks between namelist.wps and namelist.input
# for now let's take it from namelist.input and assume the two namelists are
# compatible
namelist_input_master_dict = wu.read_namelist(namelist_input_master)
max_dom = namelist_input_master_dict['&domains']['max_dom'][0]
no_of_phases = (2 + 4*(max_dom-2) + 4) + 1
phase_idx = 0
pp = 1

# control tower     
# is this a brand new run or are we continuing a previous one?
if os.path.isfile(pydown_status):
    print '\nThis is what has been already done:\n'
    phase_completed = load_status()
    ready_to_roll = False
else:
    # Initial checks to protect fragile operations
    fix_before_we_start = ''

    if os.path.isdir(wps_dir):
        wps_files = os.listdir(wps_dir)
        if not ('GRIBFILE.AAA' in wps_files and 'GRIBFILE.AAB' in wps_files):
            fix_before_we_start += "\t- There is no way for me to know which grib " \
              + "files you want linked to generate the WRF's ICs and BCs so " \
              + " please go back and link the ones you want.\n"
        if 'Vtable' not in wps_files:
            fix_before_we_start += "\t- There is no way for me to know which " \
              + "Vtable goes with your grib files of choice so please go back " \
              + "and link the one you want.\n"
    else:
        fix_before_we_start += "I couldn't find the WPS directory " \
          + wps_dir + '\n'

    if os.path.isdir(run_dir):
        run_files = os.listdir(run_dir)
        if 'namelist.input' in run_files:
            fix_before_we_start += '\t- namelist.input in your ' + run_dir \
            + ' will bother me so please move/delete it.\n'
        if 'namelist.input.master' not in run_files:
            fix_before_we_start += '\t- I expected to find a namelist.input.master ' \
              + 'in your ' + run_dir + ' please put one (possibly sensible) ' \
              + 'there for me to know what to do.\n'
        if True in ['met_em' in file for file in run_files]:
            fix_before_we_start += "\t- It looks like you've met_em files (or dir) " \
              + "in your " + run_dir + ', please move/delete them.\n'
    else:
        fix_before_we_start += "\t- I couldn't find the run directory\n" \
          + run_dir 

    if os.path.isdir(archive_dir):
        archive_files = os.listdir(archive_dir)
        if not (archive_files == ['README'] 
          or archive_files == ['README','README~']):
            # TODO this request might just be me being lazy so make it more general
            # if are so inclined
            fix_before_we_start += "\t- All I would like to see in the archive " \
              + "directory is a README file.\n" \
              + " \t\tDon't be lazy: jot down a couple of sentences in" \
              + " archive/README" \
              + " describing what you set out to achieve with this run... " \
              + " Yeah, I know it's boring, but you will thank me later ;]\n" \
              + " \t\tIf you've done this already, then get rid of any other" \
              + " file in the archive folder."
    else:
        fix_before_we_start += "\t- I couldn't find the archive directory " \
          + archive_dir + '\n'

    if fix_before_we_start != '':
        # We've got a problem
        print 'Please address the following and restart the script:'
        print fix_before_we_start
        sys.exit()
    else:
        # All is good (or the error is semantyc).
        print "Welcome to pydown. I have found all I need in order to start," \
         + "so let's begin with:"
        print phase_description(phase_idx,pp)
    
    phase_completed = [False for idx in range(no_of_phases)]
    check_with_user()
    ready_to_roll = True

if not ready_to_roll and not phase_completed[phase_idx]:
    print '\nThis is what has to be done next:\n'
    print phase_description(phase_idx,pp)
    check_with_user()   
    ready_to_roll = True
print phase_description(phase_idx,pp)
if ready_to_roll:
    os.chdir(wps_dir)
    print '\tMoved to ' + wps_dir + '.'
    # TODO I am assuming the user has already linked the correct Vtable and
    # grib files as I cannot extract that information from the namelists...

    # Let's run the wps components
    print '\tRunning geogrid.exe:'
    if 'Successful completion of geogrid' not in \
      sp.Popen('./geogrid.exe', stdout=sp.PIPE).communicate()[0]:
        # TODO handle the error in a more useful way.
        print "Houston we've got a problem... with geogrid.exe!"
    else:
        print '\t\tDone!'
    print '\tRunning ungrib.exe:'
    if 'Successful completion of ungrib' not in \
      sp.Popen('./ungrib.exe', stdout=sp.PIPE).communicate()[0]:
        # TODO handle the error in a more useful way.
        print "Houston we've got a problem... with ungrib.exe!"
    else:
        print '\t\tDone!'
    print '\tRunning metgrid.exe:'
    if 'Successful completion of metgrid' not in \
      sp.Popen('./metgrid.exe', stdout=sp.PIPE).communicate()[0]:
        # TODO handle the error in a more useful way.
        print "Houston we've got a problem... with metgrid.exe!"
    else:
        print '\t\tDone!'

    # By this stage the wps components should have executed successfully so we
    # move on to archive the appropriate metadata for future reference
    print '\tZipping metadata archive:'
    wps_metadata = 'wps_metadata.zip'
    cmd = r'zip ' + wps_metadata + ' configure.wps Vtable namelist.wps '
    # TODO I am assuming standard name for the ungrib intermediate files ->
    # this should be made more general by reading it from namelist.wps
    for file in os.listdir(os.getcwd()):
        if 'met_em' in file \
          or '.log' in file \
          or 'geo_em' in file \
          or 'GRIBFILE' in file \
          or 'FILE' in file :
            cmd += file + ' '
    # The following is to ignore both stdout and stderr from the zip command
    # TODO check if the following syntax represents a portability issue
    if sp.call(cmd.split(), stdout=open('/dev/null','w'), stderr = sp.STDOUT) == 0:
        print '\t\tDone!'

    # let's archive the metadata and move the met_em files to a met_en dir in
    # the run_dir moving them out of the wps directory will free it up for 
    # other uses
    # TODO none of the possible exceptions are caught... beware!
    print '\tArchiving metadata:'
    os.rename(wps_metadata,os.path.join(archive_dir,wps_metadata))
    print '\t\tDone!'

    print '\tMoving met_em* files to the run directory:'
    os.mkdir(os.path.join(run_dir,'met_em'))
    for file in os.listdir(os.getcwd()):
        if 'met_em' in file:
            os.rename(file,os.path.join(run_dir,'met_em',file))
    print '\t\tDone!'
    
    # let's now clean the wps_dir
    # TODO I am just being lazy here but the user 
    # should be forced to regenerate Vtable and GRIBFILE.* links every time to
    # make sure they are the right thing... beware!
          #or 'GRIBFILE' in file \
          #or 'Vtable' == file \
    print '\tCleaning the wps directory:'
    for file in os.listdir(os.getcwd()):
        if 'geo_em' in file \
          or 'FILE:' in file \
          or 'PFILE:' in file:
            os.remove(file)
    print '\t\tDone!'

    print '\tGenerating the namelist.input.d01.real:'
    current_namelist_input = \
      pdu.generate_namelist_input_d01_real(namelist_input_master,run_dir)
    print '\t\tDone!'

    print '\tLinking the just generated namelist.input:'
    os.symlink(current_namelist_input,os.path.join(run_dir,'namelist.input'))
    print '\t\tDone!'
    
    print '\tLinking the previously generated met_em files:'
    for file in os.listdir(os.path.join(run_dir,'met_em')):
        if 'd' + str(pp).zfill(2) in file:
            os.symlink(os.path.join(run_dir,'met_em',file),
              os.path.join(run_dir,file[:9] + '1' + file[10:]))
    print '\t\tDone!'

    phase_completed[phase_idx] = True
    # let's be patronizing...
    print this_is_what_you_should_do_next(phase_idx)
    check_with_user()
# we are finished let's move on to the next phase
phase_idx += 1

if not ready_to_roll and not phase_completed[phase_idx]:
    print '\nThis is what has to be done next:\n'
    print phase_description(phase_idx,pp)
    check_with_user()   
    ready_to_roll = True
print phase_description(phase_idx,pp)
if ready_to_roll:
    # let's start with a little housekeeping to clean after the real.exe
    os.chdir(run_dir)
    print '\tMoved to ' + run_dir + '.'
    print '\tMoving real log files out of the way:'
    os.mkdir(os.path.join(run_dir,'real_logs'))
    for file in os.listdir(os.getcwd()):
        if 'rsl' in file \
          or 'std' in file:
            os.rename(file,os.path.join(run_dir,'real_logs',file))
    print '\t\tDone!'
    
    print '\tGenerating the namelist.input.d01.wrf:'
    current_namelist_input = \
      pdu.generate_namelist_input_d01_wrf(namelist_input_master,run_dir)
    print '\t\tDone!'

    print '\tLinking the just generated namelist.input:'
    os.remove(os.path.join(run_dir,'namelist.input'))
    os.symlink(current_namelist_input,os.path.join(run_dir,'namelist.input'))
    print '\t\tDone!'

    phase_completed[phase_idx] = True
    print this_is_what_you_should_do_next(phase_idx)
    check_with_user()
    # we are finished let's move on to the next phase
phase_idx += 1

if not ready_to_roll and not phase_completed[phase_idx]:
    print '\nThis is what has to be done next:\n'
    print phase_description(phase_idx,pp)
    check_with_user()   
    ready_to_roll = True
print phase_description(phase_idx,pp)
if ready_to_roll:
    # let's start with a little housekeeping to clean after the wrf.exe
    os.chdir(run_dir)
    print '\tMoved to ' + run_dir + '.'
    print '\tMoving wrf log files out of the way:'
    os.mkdir(os.path.join(run_dir,'wrf_logs'))
    for file in os.listdir(os.getcwd()):
        if 'rsl' in file \
          or 'std' in file:
            os.rename(file,os.path.join(run_dir,'wrf_logs',file))
    print '\t\tDone!'

    # then move on to generating and archiving the metadata and wrfout_d01*
    print '\tZipping metadata archive:'
    wrf_metadata = 'wrf_metadata_d' + str(pp).zfill(2) + '.zip'
    cmd = r'zip -r ' + wrf_metadata + ' wrfndi_d02 ../compile.log ' \
      + '../configure.wrf wrfinput_d01 wrfbdy_d01 submit_wrf.sh ' \
      + 'submit_real.sh submit_ndown.sh namelist.input ' \
      + 'wrf_logs real_logs ndown_logs'
    # The following is to ignore both stdout and stderr from the zip command
    # TODO check if the following syntax represents a portability issue
    if sp.call(cmd.split(), stdout=open('/dev/null','w'), stderr = sp.STDOUT) == 0:
        print '\t\tDone!'

    # TODO check that the appropriate directory structure exist for the archive
    # and create it if it doesn't    
    # TODO none of the possible exceptions are caught... beware!
    print '\tArchiving metadata:'
    os.rename(wrf_metadata,os.path.join(archive_dir,wrf_metadata))
    print '\t\tDone!'

    print '\tMoving wrfout_d01* files to the archive directory:'
    for file in os.listdir(os.getcwd()):
        if 'wrfout' in file:
            os.rename(file,os.path.join(archive_dir,file))
    print '\t\tDone!'
    
    # let's now clean the wrf_dir
    print '\tCleaning the wrf directory:'
    for file_or_dir in os.listdir(os.getcwd()):
        if 'met_em.' in file_or_dir \
          or 'wrfinput' in file_or_dir \
          or 'wrfndi' in file_or_dir \
          or 'namelist.input' == file_or_dir \
          or 'wrfbdy' in file_or_dir:
            os.remove(file_or_dir)
        # only my logs directories will match 'logs'
        # TODO make this more general/robust or warn the user they should not
        # have anything that matches this check in their run_dir
        if 'logs' in file_or_dir:
            for nested_file in os.listdir(file_or_dir):
                os.remove(os.path.join(file_or_dir,nested_file))
            os.rmdir(file_or_dir)
    print '\t\tDone!'

    phase_completed[phase_idx] = True
    print this_is_what_you_should_do_next(phase_idx)
    check_with_user()
# we are finished let's move on to the next phase
phase_idx += 1

while pp < max_dom:
    pp += 1

    if not ready_to_roll and not phase_completed[phase_idx]:
        print '\nThis is what has to be done next:\n'
        print phase_description(phase_idx,pp)
        check_with_user()   
        ready_to_roll = True
    print phase_description(phase_idx,pp)
    if ready_to_roll:
        print '\tGenerating the namelist.input.d' + str(pp).zfill(2) + '.real:'
        current_namelist_input = \
          pdu.generate_namelist_input_dpp_real(pp,namelist_input_master,run_dir)
        print '\t\tDone!'

        print '\tLinking the just generated namelist.input:'
        os.symlink(current_namelist_input,os.path.join(run_dir,'namelist.input'))
        print '\t\tDone!'
        
        print '\tLinking the previously generated met_em files:'
        for file in os.listdir(os.path.join(run_dir,'met_em')):
            if 'd' + str(pp).zfill(2) in file:
                os.symlink(os.path.join(run_dir,'met_em',file),
                  os.path.join(run_dir,file[:9] + '1' + file[10:]))
        print '\t\tDone!'
        
        phase_completed[phase_idx] = True
        print this_is_what_you_should_do_next(phase_idx)
        check_with_user()
    # we are finished let's move on to the next phase
    phase_idx += 1
    
    if not ready_to_roll and not phase_completed[phase_idx]:
        print '\nThis is what has to be done next:\n'
        print phase_description(phase_idx,pp)
        check_with_user()   
        ready_to_roll = True
    print phase_description(phase_idx,pp)
    if ready_to_roll:
        # let's start with a little housekeeping to clean after the real.exe
        os.chdir(run_dir)
        print '\tMoved to ' + run_dir + '.'

        print '\tMoving real log files out of the way:'
        os.mkdir(os.path.join(run_dir,'real_logs'))
        for file in os.listdir(os.getcwd()):
            if 'rsl' in file \
              or 'std' in file:
                os.rename(file,os.path.join(run_dir,'real_logs',file))
        print '\t\tDone!'

        #let's prepare for ndown
        print '\tRename wrfinput_d01 wrfndi_02 and get rid of wrfbdy_d01:'
        os.rename('wrfinput_d01','wrfndi_d02')
        os.remove('wrfbdy_d01')
        print '\t\tDone!'

        print '\tGenerating the namelist.input.d' + str(pp).zfill(2) + '.ndown:'
        current_namelist_input = \
          pdu.generate_namelist_input_dpp_ndown(pp,namelist_input_master,run_dir)
        print '\t\tDone!'

        print '\tLinking the just generated namelist.input:'
        if os.path.isfile(os.path.join(run_dir,'namelist.input')):
            os.remove(os.path.join(run_dir,'namelist.input'))
        os.symlink(current_namelist_input,os.path.join(run_dir,'namelist.input'))
        print '\t\tDone!'

        # TODO this statement is untested
        print '\tLinking the (previously calculated) wrfout_d' \
          + str(pp-1).zfill(2)
        for file in os.listdir(archive_dir):
            if 'wrfout_d' + str(pp-1).zfill(2) in file:
                os.symlink(os.path.join(archive_dir,file),
                file[:9] + '1' + file[10:])
        print '\t\tDone!' 

        phase_completed[phase_idx] = True
        print this_is_what_you_should_do_next(phase_idx)
        check_with_user()
    # we are finished let's move on to the next phase
    phase_idx += 1
	    
    if not ready_to_roll and not phase_completed[phase_idx]:
        print '\nThis is what has to be done next:\n'
        print phase_description(phase_idx,pp)
        check_with_user()   
        ready_to_roll = True
    print phase_description(phase_idx,pp)
    if ready_to_roll:
        # let's start with a little housekeeping to clean after the ndown.exe
        os.chdir(run_dir)
        print '\tMoved to ' + run_dir + '.'

        print '\tMoving ndown log files out of the way:'
        os.mkdir(os.path.join(run_dir,'ndown_logs'))
        for file in os.listdir(os.getcwd()):
            if 'rsl' in file \
              or 'std' in file:
                os.rename(file,os.path.join(run_dir,'ndown_logs',file))
        print '\t\tDone!'

        #let's prepare for wrf
        print '\tRename wrfinput_d02 and wrfbdy_02 to wrfinput_d01 and wrfbdy_d01:'
        os.rename('wrfinput_d02','wrfinput_d01')
        os.rename('wrfbdy_d02','wrfbdy_d01')
        print '\t\tDone!'

        print '\tGenerating the namelist.input.d' + str(pp).zfill(2) + '.wrf:'
        current_namelist_input = \
          pdu.generate_namelist_input_dpp_wrf(pp,namelist_input_master,run_dir)
        print '\t\tDone!'

        print '\tLinking the just generated namelist.input:'
        os.remove(os.path.join(run_dir,'namelist.input'))
        os.symlink(current_namelist_input,os.path.join(run_dir,'namelist.input'))
        print '\t\tDone!'

        print '\tRemove the links to  wrfout_d' + str(pp-1).zfill(2) + '*'
        for file in os.listdir(run_dir):
            #if 'wrfout_d' + str(pp-1).zfill(2) in file:
            if 'wrfout_d01' in file:
                os.remove(file)
        print '\t\tDone!' 

        phase_completed[phase_idx] = True
        print this_is_what_you_should_do_next(phase_idx)
        check_with_user()
    # we are finished let's move on to the next phase
    phase_idx += 1
	    
    if not ready_to_roll and not phase_completed[phase_idx]:
        print '\nThis is what has to be done next:\n'
        print phase_description(phase_idx,pp)
        check_with_user()   
        ready_to_roll = True
    print phase_description(phase_idx,pp)
    if ready_to_roll:
        # let's start with a little housekeeping to clean after the wrf.exe
        os.chdir(run_dir)
        print '\tMoved to ' + run_dir + '.'
        print '\tMoving wrf log files out of the way:'
        os.mkdir(os.path.join(run_dir,'wrf_logs'))
        for file in os.listdir(os.getcwd()):
            if 'rsl' in file \
              or 'std' in file:
                os.rename(file,os.path.join(run_dir,'wrf_logs',file))
        print '\t\tDone!'

        # then move on to generating and archiving the metadata and wrfout_d01*
        print '\tZipping metadata archive:'
        wrf_metadata = 'wrf_metadata_d' + str(pp).zfill(2) + '.zip'
        cmd = r'zip -r ' + wrf_metadata + ' wrfndi_d02 ../compile.log ' \
          + '../configure.wrf wrfinput_d01 wrfbdy_d01 submit_wrf.sh ' \
          + 'submit_real.sh submit_ndown.sh namelist.input ' \
          + 'wrf_logs real_logs ndown_logs'
        # The following is to ignore both stdout and stderr from the zip command
        # TODO check if the following syntax represents a portability issue
        if sp.call(cmd.split(), stdout=open('/dev/null','w'), stderr = sp.STDOUT) == 0:
            print '\t\tDone!'
    
        # TODO check that the appropriate directory structure exist for the archive
        # and create it if it doesn't    
        # TODO none of the possible exceptions are caught... beware!
        print '\tArchiving metadata:'
        os.rename(wrf_metadata,os.path.join(archive_dir,wrf_metadata))
        print '\t\tDone!'
    
        print '\tMoving wrfout_d01* files to the archive directory:'
        for file in os.listdir(os.getcwd()):
            if 'wrfout' in file:
                os.rename(file,os.path.join(archive_dir,file[:8] \
                  + str(pp).zfill(2) + file[10:]))
        print '\t\tDone!'
        
        # let's now clean the wrf_dir
        print '\tCleaning the wrf directory:'
        for file_or_dir in os.listdir(os.getcwd()):
            if 'met_em.' in file_or_dir \
              or 'wrfinput' in file_or_dir \
              or 'wrfndi' in file_or_dir \
              or 'wrfbdy' in file_or_dir:
                os.remove(file_or_dir)
            if os.path.isfile(os.path.join(run_dir,'namelist.input')):
                os.remove(os.path.join(run_dir,'namelist.input'))
            # only my logs directories will match 'logs'
            # TODO make this more general/robust or warn the user they should not
            # have anything that matches this check in their run_dir
            if 'logs' in file_or_dir:
                for nested_file in os.listdir(file_or_dir):
                    os.remove(os.path.join(file_or_dir,nested_file))
                os.rmdir(file_or_dir)
        print '\t\tDone!'
    
        phase_completed[phase_idx] = True
        print this_is_what_you_should_do_next(phase_idx)
        check_with_user()
    # we are finished let's move on to the next phase
    phase_idx += 1

# the presence of a pydown.status file will be assumed to mean that we want to
# continue a partially executed run. We are now finished so we need to clean up
# the directory not to confuse the program next time it is run.
# TODO check that this exists before erasing it
os.remove(pydown_status)
print 'No more nests to run... good night and good luck.'
