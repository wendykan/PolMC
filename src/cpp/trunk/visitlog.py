# Visit 1.11.2 log file
ScriptVersion = "1.11.2"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
ShowAllWindows()
OpenDatabase("tmp.3D", 0)
# The UpdateDBPluginInfo RPC is not supported in the VisIt module so it will not be logged.
AddPlot("Volume", "0", 1, 1)
DrawPlots()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "fluenceBinNumber"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.TIFF  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK
SaveWindowAtts.width = 1024
SaveWindowAtts.height = 1024
SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
SaveWindowAtts.stereo = 0
SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
