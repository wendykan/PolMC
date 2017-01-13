"""
 Scripts to setup a multiple wavelength run on TACC ranger.
 
 $Source: /usr/data0/leipzig_work/tmat_cvs/src/polarization_MC_port/wavelength_run_setup.py,v $
"""

import os, sys, pdb
import elementtree.ElementTree as ET
from numpy import *


def wavelength_plist(dest_plist, base_plist, wavelength):
  """
  specific-wavelength plist from a base-plist
  """
  tree = ET.parse(base_plist)

  after_key = False
  for elt in tree.getiterator():
    if after_key:
      # wavelength "value" entry:
      if (elt.tag != 'real'):
        print 'elt.tag = %s' % elt.tag
        raise RuntimeError, '\"wavelength\" key with value tag not \"real\"'
      elt.text = '%f' % wavelength
      after_key = False
    elif(elt.text == 'wavelength'):
      # wavelength "key" entry:
      after_key = True

  out_fp = open(dest_plist,'wt');

  # need to _explicitely_ write the header lines to the output file (as elementtree ignores them)
  print >>out_fp, '<?xml version="1.0" encoding="UTF-8"?>'
  print >>out_fp, '<!DOCTYPE plist PUBLIC "-//Apple Computer//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">'  
  tree.write(out_fp)

  out_fp.close()
    
  
def job_script(job_filename = None, job_name = None, tasks_per_node = None, total_cores = None, 
                wc_time = None, email = None, working_directory = None, 
                NUM_THREADS = None, bin_name = None, param_filename = None):
  """
  construct a job script (as a string) from a plist-name and a run-directory name
  
  suggested input parms:
    base = <dd.mm.yy>_<run series number>
    working_directory = $WORK/run.data/<base>
    job_filename = working_directory/<base>.job
    job_name = J.<base>
    param_filename = working_directory/<plist name>
    
  remember to "mkdir -p <working_directory>"  PRIOR to running this script;
    put a copy of the plist in the working directory (use it as the <param_filename>.
    
  run this script
    
  submit the script from the working directory and that's where your log will go (which is recommended):
  cd <working directory>
  qsub <base>.job 
  """
  
  # substitution parameters 
  #   %1: <job name>  # can not start with a digit
  #   %2: <tasks/node> # integer
  #   %3: <cores total> # integer
  #   %4: <wall clock time: as <hours>:<minutes>:<seconds> >
  #   %5: <email address for notifications>
  #   %6: <working directory: preferably in the $WORK filesystem>
  #   %7: <OMP_NUM_THREADS>  # integer
  #   %8: <executable name>  # string: can include environement variables as $<environment variable>
  #   %9: <parameter file name>  # string: can include environement variables as $<environment variable>
  
  
  script_FMT = """\
#!/bin/sh

# TACC ranger submission script:

#$ -V	# Inherit the submission environment

#$ -cwd  # start job in submission directory: regardless, perhaps the "cd" will also work 

#$ -N %s  # cannot start with digit
#$ -j y	# Combine stderr and stdout
#$ -o $JOB_NAME.$JOB_ID.log	# Name of the output file

#$ -pe %dway %d	# Requests <tasks/node>way <cores total>

#$ -q normal	              # Queue name "development" (2 hours max.)
#$ -l h_rt=%s	            # Run time (hh:mm:ss) - 5 minutes
#$ -M	%s    # Use email notification address
#$ -m be	                      # Email at Begin and End of job
set -x	                        # Echo commands, use "set echo" with csh

WORKING_DIRECTORY=%s

mkdir -p $WORKING_DIRECTORY
cd $WORKING_DIRECTORY 

export OMP_NUM_THREADS=%d

# Run the MPI executable  with  parameter file:
ibrun %s %s"""
  
  out_fp = open(job_filename, 'wt')  
  print >>out_fp, script_FMT % (job_name, tasks_per_node, total_cores, wc_time, email, working_directory, NUM_THREADS, bin_name, param_filename)  
  out_fp.close()

    
def submission_script(dest_filename, run_directory_base, lambda_min, lambda_max, N_lambda):
  """
  construct a submission script (i.e. containing multiple "qsub") for a multi-wavelength set of runs
  """
  out_fp = open(dest_filename, 'wt')
  print >>out_fp, '#!/bin/sh\n'
    
  wavelength = linspace(lambda_min, lambda_max, N_lambda)
  for w in wavelength:
    print >>out_fp, 'qsub %s/lambda_%f/%s_%f.job\n' % (run_directory_base, w, run_directory_base, w)
   
  out_fp.close() 

def create_directory_structure(base_dirname, base_plist, lambda_min, lambda_max, N_lambda,
                               tasks_per_node, total_cores, wc_time, email, 
                               NUM_THREADS, bin_name):
                               
  os.makedirs(base_dirname) # error if already exists
  qsub_script_name = base_dirname + '/Q_%f_%f_%d.sh' % (lambda_min, lambda_max, N_lambda)
  qsub_fp = open(qsub_script_name,'wt'); 
  print >>qsub_fp, '#!/bin/sh\n'
                             
  for wavelength in linspace(lambda_min, lambda_max, N_lambda): 
    pos = base_dirname.rfind('/')
    if (-1 != pos):
      base = base_dirname[pos+1:]
    else:
      base = base_dirname  
    work_dir = base_dirname + '/lambda_%f' % wavelength
    plist_name = work_dir + '/' + base + '_%f.plist' % wavelength
    job_filename = work_dir + '/' + base + '_%f.job' % wavelength
    jobname = 'J.' + base + '_%f' % wavelength
    
    # create the working directory (error if already exists):
    os.makedirs(work_dir)
    # create the parameter file:
    wavelength_plist(plist_name, base_plist, wavelength)
    # create the job script:
    job_script(job_filename, jobname, 
                      tasks_per_node, total_cores, wc_time, 
                      email, work_dir, 
                      NUM_THREADS, bin_name, plist_name)
    # qsub script entry:
    print >>qsub_fp, 'qsub %s\n' % job_filename                  
  
  qsub_fp.close()                  
