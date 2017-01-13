"""
  Python command line script to merge Pol-MC output XML files to a single binary "polmc_data.world" object

  example command line
    python merge_polmc_XML.py '*.xml' world.pk.dat
"""

import os, sys, glob, pickle, pdb
import polmc_data as PD

N_required_args = 2
if (len(sys.argv) != N_required_args + 1):
  print 'usage: \"python %s \'<input glob expression>\' <output filename>\"' % sys.argv[0] 
  sys.exit(1)

glob_expr = sys.argv[1]
output_filename = sys.argv[2]

# pdb.set_trace()

print 'creating binary polmc_data: %s, from Pol-MC XML files: %s' % (output_filename, glob_expr)
world = PD.world()
for XML in glob.iglob(glob_expr):
  print '    parsing %s...' % XML
  world.read_from_XML(XML, accumulate = True)
  
fp = open(output_filename, 'wb')
pickle.dump(world, fp, protocol=pickle.HIGHEST_PROTOCOL)

fp.close()

print 'all done'
sys.exit(0)

"""
# =========== EXAMPLE: ================
# to read the polmc_data object into python:
import polmc_data as PD
import pickle

fp = open('world.pk.dat', 'rb')
world = pickle.load(fp)
fp.close()

# do something with "world"...
"""
