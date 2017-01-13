# Python script to separate photon paths into classes and save to binary (multiple runs).

import sys, os, glob, pickle, pdb
import photon_path as PP


# ============ FANCY stuff for MATPLOTLIB (aka PYLAB): =====================================================

svg_output = False

from matplotlib import rc

facecolor='(1.0,1.0,1.0)'
labelcolor='(0.0,0.0,0.0)'
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


# ================= classifying photon paths by event-type and by maximum-depth: =====================================
scattering_event_type = 1


basal_layer1 = (0.025, 0.027)
basal_layer2 = (0.025, 0.027)
basal_layer3 = (0.023, 0.027)
basal_layer4 = (0.019, 0.027)

model1_layer_depth_ranges = [(-1.0,basal_layer1[0]), (basal_layer1[0], basal_layer1[1]), (basal_layer1[1],1.0)]
model2_layer_depth_ranges = [(-1.0,basal_layer2[0]), (basal_layer2[0], basal_layer2[1]), (basal_layer2[1],1.0)]
model3_layer_depth_ranges = [(-1.0,basal_layer3[0]), (basal_layer3[0], basal_layer3[1]), (basal_layer3[1],1.0)] 
model4_layer_depth_ranges = [(-1.0,basal_layer4[0]), (basal_layer4[0], basal_layer4[1]), (basal_layer4[1],1.0)]

bin_geometry = ((-0.025, 0.025), (-0.025, 0.025), (-0.01, 0.04))
bin_shape = (200,200,200) # (Nx, Ny, Nz)


# surface basis for angle_detector:
angle_detector_perp = (sqrt(2.0)/2.0, 0.0, sqrt(2.0)/2.0)
angle_detector_para = (0.0, 1.0, 0.0)
angle_detector_normal = (-sqrt(2.0)/2.0, 0.0, sqrt(2.0)/2.0)
angle_detector_basis = (angle_detector_perp, angle_detector_para, angle_detector_normal)

# surface basis for flat_detector:
flat_detector_perp = (1.0, 0.0, 0.0)
flat_detector_para = (0.0, 1.0, 0.0)
flat_detector_normal = (0.0, 0.0, 1.0)
flat_detector_basis = (flat_detector_perp, flat_detector_para, flat_detector_normal)



UNCOLLIMATED_COSINE = 0.0
COLLIMATED_COSINE = 0.95
N_PHOTONS = 1.28e9



def output_detector_type(base_directory, MODEL_LAYER_DEPTHS, DETECTOR_TYPE, DETECTOR_BASIS, UNCOLLIMATED_COSINE, COLLIMATED_COSINE, N_PHOTONS):

  figure_directory = base_directory + '/fig_tmp'
  if not os.path.exists(figure_directory):
    os.makedirs(figure_directory) # error if already exists...

  pos = base_directory.rfind('/')
  if (pos != -1):
    leaf = base_directory[pos+1:]
  else:
    leaf = base_directory  

  photons_filename = base_directory + '/' + DETECTOR_TYPE + '_detector_total_photons.dat'

  print 'separating photon-paths into classes (directory: %s, detector: %s)...' % (leaf, DETECTOR_TYPE)
  photon_classes = PP.photon_depth_classes(PP.photon_paths_iter(photons_filename), MODEL_LAYER_DEPTHS, scattering_event_type)

  for n_class, classified_photons in enumerate(photon_classes):
    print 'generating fluence bins for model %s; detector: %s, class: %d' % (leaf, DETECTOR_TYPE, n_class)

    fluence3D, detector_fluence = PP.volume_bins(classified_photons, 
                                       bin_geometry, bin_shape, 
                                       'fluence', False, 1.0, 0.0, 1,
                                       DETECTOR_BASIS[0], DETECTOR_BASIS[1], DETECTOR_BASIS[2], 
                                       UNCOLLIMATED_COSINE, maxAcceptanceCosine = 1.0)
    output_filename_base = figure_directory + '/' + DETECTOR_TYPE + '_class_' + str(n_class) + '_uncollimated_fluence'
    fp = open(output_filename_base + '.pk.dat', 'wb')
    pickle.dump({'fluence3D':fluence3D, 'detector_fluence':detector_fluence, 'N_photons':N_PHOTONS}, fp)
    fp.close

    imshow((fluence3D[:, fluence3D.shape[1]/2, :]).T, aspect='equal'); colorbar()
    savefig(output_filename_base + '.png', format='png', dpi=300)
    close()                                   


    fluence3D, detector_fluence = PP.volume_bins(classified_photons, 
                                       bin_geometry, bin_shape, 
                                       'fluence', False, 1.0, 0.0, 1,
                                       DETECTOR_BASIS[0], DETECTOR_BASIS[1], DETECTOR_BASIS[2], 
                                       COLLIMATED_COSINE, maxAcceptanceCosine = 1.0)
    output_filename_base = figure_directory + '/' + DETECTOR_TYPE + '_class_' + str(n_class) + '_collimated_fluence'
    fp = open(output_filename_base + '.pk.dat', 'wb')
    pickle.dump({'fluence3D':fluence3D, 'detector_fluence':detector_fluence, 'N_photons':N_PHOTONS}, fp)
    fp.close

    imshow((fluence3D[:, fluence3D.shape[1]/2, :]).T, aspect='equal'); colorbar()
    savefig(output_filename_base + '.png', format='png', dpi=300)
    close()                                   
  
# ---------- model1: -----------------------------------------------------------------
base_directory = '/mnt/hd/POLMC.07.2010/run.data/28.07.2010_1' 
MODEL_LAYER_DEPTHS = model1_layer_depth_ranges

# output_detector_type(base_directory, MODEL_LAYER_DEPTHS, 'angle', angle_detector_basis, UNCOLLIMATED_COSINE, COLLIMATED_COSINE, N_PHOTONS)
# output_detector_type(base_directory, MODEL_LAYER_DEPTHS, 'flat', flat_detector_basis, UNCOLLIMATED_COSINE, COLLIMATED_COSINE, N_PHOTONS)

# ---------- model2: -----------------------------------------------------------------
base_directory = '/mnt/hd/POLMC.07.2010/run.data/28.07.2010_2' 
MODEL_LAYER_DEPTHS = model2_layer_depth_ranges

output_detector_type(base_directory, MODEL_LAYER_DEPTHS, 'angle', angle_detector_basis, UNCOLLIMATED_COSINE, COLLIMATED_COSINE, N_PHOTONS)
# output_detector_type(base_directory, MODEL_LAYER_DEPTHS, 'flat', flat_detector_basis, UNCOLLIMATED_COSINE, COLLIMATED_COSINE, N_PHOTONS)

# ---------- model3: -----------------------------------------------------------------
base_directory = '/mnt/hd/POLMC.07.2010/run.data/28.07.2010_3' 
MODEL_LAYER_DEPTHS = model3_layer_depth_ranges

output_detector_type(base_directory, MODEL_LAYER_DEPTHS, 'angle', angle_detector_basis, UNCOLLIMATED_COSINE, COLLIMATED_COSINE, N_PHOTONS)
# output_detector_type(base_directory, MODEL_LAYER_DEPTHS, 'flat', flat_detector_basis, UNCOLLIMATED_COSINE, COLLIMATED_COSINE, N_PHOTONS)

# ---------- model4: -----------------------------------------------------------------
base_directory = '/mnt/hd/POLMC.07.2010/run.data/28.07.2010_4' 
MODEL_LAYER_DEPTHS = model4_layer_depth_ranges

output_detector_type(base_directory, MODEL_LAYER_DEPTHS, 'angle', angle_detector_basis, UNCOLLIMATED_COSINE, COLLIMATED_COSINE, N_PHOTONS)
# output_detector_type(base_directory, MODEL_LAYER_DEPTHS, 'flat', flat_detector_basis, UNCOLLIMATED_COSINE, COLLIMATED_COSINE, N_PHOTONS)
