import os
import sys
import cPickle

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
              + 'We start with running wps and doing the housekeeping ' \
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
        
def check_with_user():
    answer = raw_input('Whaddoyouwonnadonow?[c,s]:')
    if answer not in ['c','s']:
        check_with_user()
    if answer == 's':
        dump_status()
        sys.exit('Status saved: see yaa!')

def dump_status():
    cPickle.dump((phase_idx,pp),open(pydown_status,'w'))
    return

def load_status():
    return cPickle.load(open(pydown_status,'r'))
    


#TODO this will be read from one of the namelists
max_dom = 3
# TODO choose what to assume about where pydown will be run from
# in the first instance
#pydown_dir = os.getcwd()
pydown_dir = '/nfs/1/home/vbisign/wrf/pydown_sbox'
pydown_status = os.path.join(pydown_dir,'pydown.status')
wps_dir = os.path.join(pydown_dir,'WPS')
run_dir = os.path.join(pydown_dir,'WRFV2','run')

# control tower     
# is this a brand new run or are we continuing a previous one?
if os.path.isfile(pydown_status):
    print 'This is what has been already done:'
    old_phase_idx, old_pp = load_status()
    ready_to_roll = False
else:
    print 'Brand new'
    ready_to_roll = True
phase_idx = 0
pp = 1


if not ready_to_roll and phase_idx == old_phase_idx:
    ready_to_roll = True
print phase_description(phase_idx,pp)
if ready_to_roll:
    check_with_user()
# we are finished let's move on to the next phase
phase_idx += 1

if not ready_to_roll and phase_idx == old_phase_idx:
    ready_to_roll = True
print phase_description(phase_idx,pp)
if ready_to_roll:
    check_with_user()
    # we are finished let's move on to the next phase
phase_idx += 1

if not ready_to_roll and phase_idx == old_phase_idx:
    ready_to_roll = True
print phase_description(phase_idx,pp)
if ready_to_roll:
    check_with_user()
# we are finished let's move on to the next phase
phase_idx += 1

while pp < max_dom:
    pp += 1

    if not ready_to_roll and phase_idx == old_phase_idx:
        ready_to_roll = True
    print phase_description(phase_idx,pp)
    if ready_to_roll:
        check_with_user()
    # we are finished let's move on to the next phase
    phase_idx += 1
    
    if not ready_to_roll and phase_idx == old_phase_idx:
        ready_to_roll = True
    print phase_description(phase_idx,pp)
    if ready_to_roll:
        check_with_user()
    # we are finished let's move on to the next phase
    phase_idx += 1
	    
    if not ready_to_roll and phase_idx == old_phase_idx:
        ready_to_roll = True
    print phase_description(phase_idx,pp)
    if ready_to_roll:
        check_with_user()
    # we are finished let's move on to the next phase
    phase_idx += 1
	    
    if not ready_to_roll and phase_idx == old_phase_idx:
        ready_to_roll = True
    print phase_description(phase_idx,pp)
    if ready_to_roll:
        check_with_user()
    # we are finished let's move on to the next phase
    phase_idx += 1

# the presence of a pydown.status file will be assumed to mean that we want to
# continue a partially executed run. We are now finished so we need to clean up
# the directory not to confuse the program next time it is run.
os.remove(pydown_status)
print 'No more nests to run... good night and good luck.'
