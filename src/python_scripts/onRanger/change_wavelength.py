"""
  Change all "wavelength" key-value pair elements in a plist to match the wavelength specified on the command line.
  
  example command line:
    python change_wavelength.py <input plist> <wavelength> <output plist>
  
  $Source: /usr/data0/leipzig_work/tmat_cvs/src/polarization_MC_port/change_wavelength.py,v $
"""

import os, sys, pdb
import elementtree.ElementTree as ET

if (len(sys.argv) != 4):
  print 'usage: python change_wavelength.py <input file> <wavelength> <output file>'
  sys.exit(1)

try:

  input_filename =  sys.argv[1]
  wavelength = float(sys.argv[2])
  output_filename = sys.argv[3]

  tree = ET.parse(input_filename)

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
  
  # need to _explicitely_ write the header lines to the output file (as elementtree ignores them)
  out_fp = open(output_filename,'wt');

  print >>out_fp, '<?xml version="1.0" encoding="UTF-8"?>'
  print >>out_fp, '<!DOCTYPE plist PUBLIC "-//Apple Computer//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">'  
  tree.write(out_fp)
  out_fp.close()
  
except:
  print 'change_wavelength.py ERROR parsing plist file %s to new file %s' % (input_filename, output_filename)
  sys.exit(1)
  
sys.exit(0)  
