# This script lets you choose which attribute you want to plot in 3D and saves a TIF image
# To run this script do: > visit -cli -s wendyPlotScript.py
# Wendy Kan 7/18/10

import os,sys

filetype = int(raw_input( "which 3D plot do you want to plot? 1: energyBinNumber, 2: fluenceBinNumber "))

if filetype == 1:
	arg_name = 'energyBinNumber'
else:
	arg_name = 'fluenceBinNumber'

os.system( 'xpath out_0.xml "/simulation/' + str(arg_name) + '" > tmp.3D')


OpenDatabase('tmp.3D',0) 
AddPlot("Volume", "0", 1, 1)
DrawPlots()

s = SaveWindowAttributes() 
s.format = s.TIFF
s.fileName = arg_name 
#s.width, s.height = 1024,768 
#s.screenCapture = 0 
SetSaveWindowAttributes(s) 

SaveWindow()

sys.exit()