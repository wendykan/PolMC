"""
  Calculation of detector position for multilayer structure.
  
  $Source: $
"""

from numpy import *

# Starting with the layer containing the detector, specify layer structure as a list of (<thickness>, <refractive index>) tuples:

# sourceLayer, epithelium, basal, stromal_membrane, stroma
layers = [(1.0, 1.335), (0.025, 1.369), (0.002, 1.369), (0.003, 1.41), (0.1,1.369)]

def transverse_offset(angle, z_offset, layers):
  """
  With specified source angle and layers structure, calculate the transverse offset (e.g. the offset in the xy-plane),
  corresponding to a ray propagation through the layers stack.
  
  input parameters:
    angle: angle of detector optical axis with respect to the vertical axis
    z_offset: vertical offset of detector from first layer interface
    layers:  list of (<thickness>, <refractive index>) tuples starting from the layer containing the detector
    
  returns:
    transverse offset   
  """
  
  # initial detector offset
  x_offset = z_offset * sin(angle)
  
  sin_i = sin(angle)
  n_i = layers[0][1]
  # repeatedly apply Snell's law: n_i sin(\theta_i) = n_t sin(\theta_t):
  for layer in layers[1:]:
    z_t = layer[0]
    n_t = layer[1]
    sin_t = n_i*sin_i/n_t
    theta_t = arcsin(sin_t)
    x_offset += z_t * tan(theta_t)
    # prepare for next layer
    sin_i = sin_t
    n_i = n_t
    
  return x_offset
  
