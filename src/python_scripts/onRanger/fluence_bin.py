"""
python script to extract fluence bin information from Pol-MC xml output:

  INSTALL the following: NUMPY, SCIPY, MATPLOTLIB, ELEMENTTREE
"""
output_file = 'out_0.xml'

import elementtree.ElementTree as ET
import StringIO as STR
from numpy import *
from pylab import *
import sys

input_file = sys.argv[1]
tree = ET.parse(sys.argv[1])

if (len(sys.argv) > 2):
  output_file = sys.argv[2]
else:
  output_file = input_file[0:(input_file).rfind('.')] + ".png"

if (len(sys.argv) > 3):
  if (sys.argv[3] == 'log'):
    bLog = True
else:
  bLog = False

# string containing lines: nx ny nz <value>
fluenceBins = tree.find('/fluenceBinNumber')
Nx = float(fluenceBins.attrib['Nx'])
Ny = float(fluenceBins.attrib['Ny'])
Nz = float(fluenceBins.attrib['Nz'])

istr = STR.StringIO(fluenceBins.text)
bins = []
for line in istr:
  tokens = line.split()
  if len(tokens) >= 4:
    bins.append(float(tokens[3]))
istr.close()

aBins = array(bins).reshape((Nx,Ny,Nz))

if not bLog:
  imshow(aBins[:,Ny/2.0,:]); colorbar(); 
else:
  imshow(log10(aBins[:,Ny/2.0,:]+finfo(double).eps)); colorbar();  
 
savefig(output_file, dpi=300, format="png")
