"""
  script to test "polmc_data" module
"""

import sys, os, glob, pickle, pdb

import polmc_data as PD

base_directory = '/usr/data0/leipzig_work/working_copy/tmp/POLMC_BUILD/test/26.07.2010_Dy_1'

world = PD.world()
for n, XML in enumerate(glob.glob(base_directory + '/*.xml')):
  print 'accumulating from XML file number %d (%s)' % (n, XML)
  world.read_from_XML(XML,accumulate = True)

pdb.set_trace()
# PRINT-OUT a the entire data object (this is an easy way to see what all of the key names are, and how the dict are nested):  
fp = open('world_print.txt', 'wt')
print >>fp, '======================================== WORLD: ================================================='
print >>fp, world
print >>fp, world['objects']
print >>fp, '================================================================================================='
fp.close()

# SAVE the data object to a file (this is really useful, so you don't have to re-accumulate it again):
fp = open('world.pk.dat','wb')
pickle.dump(world, fp, protocol=pickle.HIGHEST_PROTOCOL)
fp.close()

# YOU CAN LOAD the data object again like this:
#
# import pickle
# fp = open('world.pk.dat','rb')
# world = pickle.load(fp)
# fp.close()
#
