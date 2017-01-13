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
imshow(aBins[:,Ny/2.0,:]); colorbar(); 
 
savefig(output_file, dpi=300, format="png")
