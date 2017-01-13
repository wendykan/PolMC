"""
  Implement a particle size (or other property) mixture model for polarization Monte Carlo (Pol-MC) using multiple input "plist" files.

  $Source: /usr/data0/leipzig_work/tmat_cvs/src/polarization_MC_port/mixture_model.py,v $
"""
import scipy
from scipy import stats


def split_plists(dest_dir, src_name, key, dist, dist_parms):
  # create directory of mixture-model plists from source plist
  
def merge_output_xml(dest_name, key_list, src_dir):
  # create merged (i.e. accumulated) output XML file from a directory containing multiple output XML files  

def split(dest_dir, src_name, key, dist, dist_parms):
  """
  split a source plist into a mixture-model mapped set of destination plists:
    input parms:
      dest_dir: destination directory
      src_name: filename of source plist
      key: key name for parameter to map to mixture model (e.g. <particle radius>)
      dist: distribution name to use for the mapping
      dist_parms: additional distribution parameters (as a list or tuple)
        (first parameter will be the value of the parameter being mapped as the mean, and initially this should just be a "placeholder" parm)        
  """

  % read the source XML file:

  % generate N interval tuples over the real line, each with equal probability, 
  %   and calculate the associated incomplete expectations:
  dist_parms[0] = src[key]
  intervals, values = balance_partial_cdf(dist, N, dist_parms)

  for n in range(N):
    % initialize destination as a copy:
    dest = src.copy()
    
    % modify the parameter:
    dest[key] = values[n]
  
    % generate destination plist-name (as fully qualified path):
    dest_name = '%s/%s_%d.xml' % (dest_dir, src, nplist)  % ACTUALLY generate as "base_00nnn.xml" with constant width and leading zeros

    % output the destination plist as XML:
  
  
def merge(dest, key_list, src):
  """
  accumulate the values associated with the key list from the source to the destination XML
    input parms:
      dest: destination dictionary
      key_list: list of keys to accumulate
      src: source dictionary
  """
  if (empty(dest)):
    dest = src
  else:
    for key in key_list:
      dest[key] += src[key]


def intervals, values = balance_partial_cdf(dist, N, dist_parms):
  """
  generate N interval tuples over the real line, each with equal probability, 
     and calculate the associated incomplete expectations:
    input parms:
      dist: name of distribution to use for the mapping
      N: number of sub-intervals
      dist_parms: parameters for the distribution
  """
  
