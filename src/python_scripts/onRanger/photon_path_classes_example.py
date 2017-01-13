"""
  Examples for working with binary photon-paths files:
  
    -- classify photon-paths: create lists for each photon-path class;
    
    -- display a fluence image from a list (or other iterable) of photon-paths,
         optionally restrict acceptance cosine at the detector;
"""

import sys, os, glob, pickle, pdb
import photon_path as PP


# ============ FANCY stuff for MATPLOTLIB (aka PYLAB): =====================================================

svg_output = False

from matplotlib import rc

facecolor='(0.0,0.0,0.2)'
labelcolor='(0.8,0.8,1.0)'
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
base_directory = '/usr/data1/polmc_project/run.data/28.07.2010_3'  # set this to YOUR base directory
figure_save_directory = '/usr/data1/polmc_project/fig_tmp'

# forward surface normal for angle_detector:
detector_normal = (-sqrt(2.0)/2.0, 0.0, sqrt(2.0)/2.0)

test_filename = base_directory + '/angle_detector_total_photons.dat'
# ********** THIS IS MODEL3 ***********

# =========================================================================================================

# ================= classifying photon paths by event-type and by maximum-depth: =====================================
scattering_event_type = 1

# -- the classifier function actually _additionally_ selects by event-type:	--
basal_layer1 = (0.025, 0.027)
basal_layer2 = (0.025, 0.027)
basal_layer3 = (0.023, 0.027)  # FOR THIS EXAMPLE, we use THIS ONE.
basal_layer4 = (0.019, 0.027)

model1_layer_depth_ranges = [(-1.0,basal_layer1[0]), (basal_layer1[0], basal_layer1[1]), (basal_layer1[1],1.0)]
model2_layer_depth_ranges = [(-1.0,basal_layer2[0]), (basal_layer2[0], basal_layer2[1]), (basal_layer2[1],1.0)]
model3_layer_depth_ranges = [(-1.0,basal_layer3[0]), (basal_layer3[0], basal_layer3[1]), (basal_layer3[1],1.0)] # AND THIS one
model4_layer_depth_ranges = [(-1.0,basal_layer4[0]), (basal_layer4[0], basal_layer4[1]), (basal_layer4[1],1.0)]

bin_geometry = ((-0.025, 0.025), (-0.025, 0.025), (-0.01, 0.04))
bin_shape = (200,200,200) # (Nx, Ny, Nz)


# surface basis for angle_detector:
angle_detector_perp = (sqrt(2.0)/2.0, 0.0, sqrt(2.0)/2.0)
angle_detector_para = (0.0, 1.0, 0.0)
angle_detector_normal = (-sqrt(2.0)/2.0, 0.0, sqrt(2.0)/2.0)


# surface basis for flat_detector:
flat_detector_perp = (1.0, 0.0, 0.0)
flat_detector_para = (0.0, 1.0, 0.0)
flat_detector_normal = (0.0, 0.0, 1.0)


# classify photon-paths for model3:
print 'separating photon-paths into classes:'
photon_classes = PP.photon_depth_classes(PP.photon_paths_iter(test_filename), model3_layer_depth_ranges, scattering_event_type)

# bin volume fluence for photon-paths with maximum depth inside the basal layer:
print 'calculating volume fluence (detector surface normal: %s, acceptance-cosine in (%f, %f)):' %\
         
volume_fluence, detector_fluence = PP.volume_bins(photon_classes[1], 
                                     bin_geometry, bin_shape, 
                                     'fluence', False, 1.0, 0.0, 1,
                                     angle_detector_perp, angle_detector_para, angle_detector_normal, 
                                     minAcceptanceCosine = 0.0, maxAcceptanceCosine = 1.0)

# --------------- display the xz-plane: -----------------------------------------------
imshow(volume_fluence[:, volume_fluence.shape[1]/2, :]); colorbar()
title('fluence (class I: basal layer photons)')
savefig(figure_save_directory + '/basal_layer_fluence' + '.png', format='png', dpi=300); 
if svg_output:
  savefig(figure_save_directory + '/basal_layer_fluence' + '.svg', format='svg');
close()     
                                       

# -------------- recalculate the volume fluence using a narrower acceptance angle: ------------

# bin volume fluence for photon-paths with maximum depth inside the basal layer:
print 'calculating volume fluence (detector surface normal: %s, acceptance-cosine in (%f, %f)):' %\
         (detector_normal, 0.95, 1.0)
volume_fluence, detector_fluence = PP.volume_bins(photon_classes[1], 
                                     bin_geometry, bin_shape, 
                                     'fluence', False, 1.0, 0.0, 1,
                                     angle_detector_perp, angle_detector_para, angle_detector_normal, 
                                     minAcceptanceCosine = 0.95, maxAcceptanceCosine = 1.0)


# --------------- display the xz-plane: -----------------------------------------------
imshow(volume_fluence[:, volume_fluence.shape[1]/2, :]); colorbar()
title('collimated-detector fluence (class I: basal layer photons)')
savefig(figure_save_directory + '/basal_layer_collimated_fluence' + '.png', format='png', dpi=300); 
if svg_output:
  savefig(figure_save_directory + '/basal_layer_collimated_fluence' + '.svg', format='svg');
close()     
                                       

