"""
  Create volume fluence image from "polmc_data"object.
  
"""

import sys, os, glob, pickle, pdb

# ============ FANCY stuff for MATPLOTLIB (aka PYLAB): =====================================================

svg_output = False

from matplotlib import rc

#facecolor='(0.0,0.0,0.2)'
#labelcolor='(0.8,0.8,1.0)'
facecolor='(1.0, 1.0, 1.0)'
labelcolor='(0.0, 0.0, 0.0)'

if (not svg_output):
  # rc("font",**{"family":"sans-serif","sans-serif":["Helvetica"]})
  rc("text", usetex=True)
rc("font",size=14)
rc("font",size=14)
rc("figure",facecolor=facecolor)
rc("savefig",facecolor=facecolor)
rc("axes",labelcolor=labelcolor)
rc("text",color=labelcolor)
rc("xtick",color=labelcolor)
rc("ytick",color=labelcolor)

from numpy import *; from pylab import *

# ==========================================================================================================


# ============================== SET these values to what you NEED: =======================================
base_directory = '/Volumes/ELEMENTS/wendykan/polmc_project/run.data/28.07.2010_1'  # set this to YOUR base directory
figure_save_directory = '/Users/wendykan/Dropbox/polmc_paper2010/polmc_draft2/figures'
minAcceptanceCosine = 0.0
# =========================================================================================================

# this just automatically generates the name of the file for the pickled "polmc_data" object from the base
#   directory name (i.e. all of it is optional):
pos = base_directory.rfind('/')
if (pos != -1):
  leaf = base_directory[pos+1:]

  # load the pickled polmc_data object:  
  fp = open(base_directory + '/world.' + leaf + '.pk.dat')
  world = pickle.load(fp)
  fp.close()

  volume_fluence = world['fluenceBinNumber']
else:
  raise RuntimeError, 'total_fluence_example.py: unable to determine file name for binary data'


nm = 1.0e-7 # cm

# --------------- display the xz-plane: -----------------------------------------------
A = volume_fluence[:, volume_fluence.shape[1]/2, :]
A = A.T
imshow(A); colorbar()
title('total volume fluence')
savefig(figure_save_directory + '/total_fluence' + '.png', format='png', dpi=300); 
if svg_output:
  savefig(figure_save_directory + '/total_fluence' + '.svg', format='svg');
close()     
                                       
